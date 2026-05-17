"""
de_quackling: Database manager for differential expression analysis data.

Small wrapper around a DuckDB file used to store and query experimental and
gene result data. Key behaviors:

- Creates and manages tables (`experimental_data`, `gene_results`, `genes`).
- Provides helper methods to ingest and normalize incoming dataframes.
- Query and delete helpers accept simple filter kwargs like `pvalue__gt=0.05`.
- If a filter column doesn't match a real column exactly, a close-match
    suggestion is raised (via difflib). If there is no close match but the
    `other_info` JSON column exists, the filter is applied against that JSON
    using `other_info ->> '<key>'`. Numeric comparisons use `CAST(... AS DOUBLE)`.

Workflow summary (high level):
1. Open connection with `with DBManager() as db:` which creates tables/indexes.
2. Use `insert_to_database()` to normalize and insert a DataFrame (preprocess_df).
3. `query()` and delete methods accept keyword filters, resolve them to
     concrete SQL clauses, and execute against DuckDB.
"""
import duckdb
import re
import pandas as pd
import os
import json
from pathlib import Path
import hashlib
import sys
import difflib
from .exceptions import ProcessingError, DuplicateExperimentError, DuplicateGeneTableError
from .exceptions import ProcessingError, DuplicateExperimentError, DuplicateGeneTableError

class de_quackling:
    def __init__(self, db_path='SQL.duckdb'):
        self.db_path = db_path
        self.conn = None
        self.insertion_columns = {
            'experiment_id': ['experiment_id', 'exp_id', 'id'],
            'gene_symbol':['symbol', 'gene_symbol'],
            'ensembl_id':['ensembl_id', 'ensemblid', 'ensembl_gene_id', 'gene_id'],
            'gene_name': ['gene_name', 'gene', 'genename'],
            'log2fc': ['log2fc', 'log2foldchange', 'log2fold'],
            'logCPM': ['logcpm', 'basemean', 'logcpm'],
            'pvalue': ['pvalue', 'p-value', 'p_value'],
            'padj': ['padj', 'fdr', 'false_discovery_rate'],
            'other_info': ['other_info', 'extra_info', 'json_info'],
            
        }
        
        
    def __enter__(self):
        self.conn = duckdb.connect(self.db_path)

        # Maintain a sequence for experiment primary keys.
        self.conn.execute('''CREATE SEQUENCE IF NOT EXISTS seq_experiment_id START 1''')

        # Create tables if they don't already exist.
        self.conn.execute('''CREATE TABLE IF NOT EXISTS experimental_data
                            (experiment_id INTEGER PRIMARY KEY DEFAULT nextval('seq_experiment_id'), 
                             tool VARCHAR, date VARCHAR, file VARCHAR, 
                             experiment_name VARCHAR, comparison_label VARCHAR, data_signature VARCHAR)''')
        
        
        self.conn.execute('''CREATE SEQUENCE IF NOT EXISTS seq_gene_id START 1''')
        
        self.conn.execute('''CREATE TABLE IF NOT EXISTS gene_results
                            (id INTEGER PRIMARY KEY DEFAULT nextval('seq_gene_id'),
                            experiment_id INTEGER,
                            gene_symbol VARCHAR,
                            ensembl_id VARCHAR,
                            log2fc DOUBLE,
                            logCPM DOUBLE,
                            pvalue DOUBLE,
                            padj DOUBLE,
                            other_info VARCHAR,
                            FOREIGN KEY (experiment_id) REFERENCES experimental_data(experiment_id))
                            ''')
        
        self.conn.execute('''CREATE TABLE IF NOT EXISTS genes
                            (
                            symbol VARCHAR, 
                            id VARCHAR,
                            ensembl_id VARCHAR,
                            alias_symbol VARCHAR[],
                            prev_symbol VARCHAR[], species VARCHAR,
                            PRIMARY KEY (id, species) 

                            )''')
        self.conn.execute("CREATE INDEX IF NOT EXISTS idx_ref_symbol ON genes (symbol)")
        self.conn.execute("CREATE INDEX IF NOT EXISTS index_ref_ensembl_id ON genes (ensembl_id)")

        
        
        
        # Create indexes
        self.conn.execute('CREATE INDEX IF NOT EXISTS idx_gene_name_lookup ON gene_results(gene_symbol, experiment_id)')
        
        # Get table and column information
        self.tables = self.conn.execute("SELECT table_name FROM information_schema.tables WHERE table_schema='main'").fetchall()
        self.gene_columns = self.conn.execute("SELECT * FROM information_schema.columns WHERE table_name='gene_results'").fetchall()
        self.experiment_columns = self.conn.execute("SELECT * FROM information_schema.columns WHERE table_name='experimental_data'").fetchall()
        
        
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        if self.conn:
            self.conn.close()
    
    def connect(self):
        if not self.conn:
            self.__enter__()
        return self

    def _split_filter_key(self, filter_key):
        # Parse the optional operator suffix from filter keys like `pvalue__gt`.
        if '__' in filter_key:
            filter_column, _, suffix = filter_key.rpartition('__')
            if suffix in ['gt', 'lt', 'gte', 'lte', 'ne']:
                return filter_column, suffix
        return filter_key, None

    def _resolve_filter_column(self, filter_column, actual_columns, table):
        # if the filter column is in the actual columns, returns the column as is
        lower_columns = [col.lower() for col in actual_columns]
        if filter_column.lower() in lower_columns:
            return actual_columns[lower_columns.index(filter_column.lower())]

        # If the user typoed a column, proactively suggest the closest real name.
        close_matches = difflib.get_close_matches(filter_column.lower(), lower_columns, n=1, cutoff=0.8)
        if close_matches:
            suggested = actual_columns[lower_columns.index(close_matches[0])]
            raise ProcessingError(f"Column '{filter_column}' does not exist in table '{table}'. Did you mean '{suggested}'?")

        # Allow fallback to other_info JSON if the requested field is not a native column.
        if 'other_info' in actual_columns:
            return None

        raise ProcessingError(
            f"Column '{filter_column}' does not exist in table '{table}'. Available columns: {', '.join(actual_columns)}"
        )

    def _build_filter_clause(self, filter_column, operator, value, actual_columns, table):
        """Return an SQL clause and params for the given filter.

        If the filter_column maps to a real column, use that column for the
        comparison. If it does not but `other_info` exists, translate to a
        JSON extraction using the `->>` operator. Numeric comparisons are
        supported by casting the extracted JSON text to DOUBLE.
        """
        actual_column = self._resolve_filter_column(filter_column, actual_columns, table)
        operators = {
            'gt': '>',
            'lt': '<',
            'gte': '>=',
            'lte': '<=',
            'ne': '!='
        }

        if actual_column is not None:
            if actual_column in ['experiment_name', 'comparison_label', 'gene_symbol', 'tool']:
                return f'{actual_column} ILIKE ?', f'%{value}%'
            if operator is None:
                return f"{actual_column} = ?", value
            if actual_column=='date':
                value=str(pd.to_datetime(value).strftime('%Y-%m-%d'))

            if operator not in operators:
                raise ProcessingError(f"Invalid operator '{operator}' for column '{filter_column}'.")

            return f"{actual_column} {operators[operator]} ?", value

        if 'other_info' not in actual_columns:
            raise ProcessingError(
                f"Column '{filter_column}' does not exist in table '{table}', and no JSON fallback is available."
            )

        # Use ->> to extract the top-level JSON field as text from other_info.
        print(f'{filter_column} does not exist in {table}, querying JSON...')
        json_key = filter_column.replace("'", "''")
        json_expr = f"other_info ->> '{json_key}'"

        if operator is None:
            return f"({json_expr}) = ?", value
        if operator == 'ne':
            return f"({json_expr}) != ?", value

        # For numeric comparisons, CAST the extracted text to DOUBLE
        if operator in ('gt', 'lt', 'gte', 'lte'):
            return f"CAST({json_expr} AS DOUBLE) {operators[operator]} ?", value

        raise ProcessingError(
            f"JSON filtering only supports equality, inequality, and numeric comparisons for '{filter_column}'."
        )

    def _find_date_and_file(self, info):
        if os.path.isfile(os.path.abspath(info)):
            date_pattern = re.compile(r'(\d{4}-\d{2}-\d{2})') 
            match = date_pattern.search(info) # first tries to find the date using regex on the file name
            file_path=os.path.abspath(info)
            if match:
                return pd.to_datetime(match.group(1), format='%Y-%m-%d'), file_path
            else: # if the search doesn't work, then uses the path library to find the time, or else uses none
                path=Path(file_path)
                if path.stat().st_mtime:
                    return pd.to_datetime(path.stat().st_mtime, unit='s').strftime('%Y-%m-%d'), file_path
        else:
            file_path='dataframe'
            return pd.Timestamp.now().strftime('%Y-%m-%d'), file_path
    
    def _get_data_signature(self, df):
        if isinstance(df, pd.DataFrame):
            pass
        elif isinstance(df, list) and isinstance(df[0], dict):
            df=pd.DataFrame(df)
        elif os.path.exists(os.path.abspath(df)):
            df=self._preprocess_df(df, 0)
        else:
            raise ProcessingError(f'df must be a pandas DataFrame, a list of dictionaries, or a file path, not {type(df)}')
        sample_data = df[['gene_symbol', 'pvalue', 'log2fc']].head(100).to_string()
    
        # Generate a unique string (hash) from that data
        return hashlib.md5(sample_data.encode()).hexdigest()

    def _preprocess_df(self, info, id):
        """Normalize incoming data into a dataframe suitable for insertion.

        Accepts a DataFrame, list-of-dicts, or a filepath. Normalization steps:
        - Convert supported inputs to a pandas DataFrame
        - Rename known input column variants to canonical insertion column names
        - Ensure expected columns exist and aggregate any extra columns into
          the `other_info` JSON column (as dicts)
        Returns a DataFrame with columns: experiment_id, gene_name, log2fc,
        logCPM, pvalue, padj, other_info
        """

        if isinstance(info, list):
            if len(info) == 0:
                return
            if isinstance(info[0], dict):
                info = pd.DataFrame(info)
            elif isinstance(info[0], pd.Series):
                info = pd.DataFrame(info)
            else:
                raise ProcessingError(f"Info list contains unsupported element type {type(info[0])}.")
        elif isinstance(info, pd.DataFrame):
            info = info.copy()
        elif os.path.isfile(os.path.abspath(info)):
            info = pd.read_csv(os.path.abspath(info), sep=None, engine='python')
        else:
            raise ProcessingError(f'Info is in an unexpected format {type(info)}; please provide a pandas DataFrame, a list of dictionaries, or a file path to a CSV or TSV file.')
        
        # Map common source column variants to canonical insert columns.
        rename_map = {}
        flat_map = {v.lower(): k for k, variants in self.insertion_columns.items() for v in variants}
        rename_dict = {col: flat_map[col.lower()] for col in info.columns if col.lower() in flat_map}
        info.rename(columns=rename_dict, inplace=True)

        if 'ensembl_id' in info.columns and 'gene_symbol' not in info.columns:
            info['gene_symbol'] = None
        if 'gene_symbol' in info.columns and 'ensembl_id' not in info.columns:
            info['ensembl_id'] = None


        if 'gene_name' in info.columns:
            if info['ensembl_id'].notnull().any(): # ensembl_id already exists
                info['gene_symbol'] = info['gene_name']
            elif info['gene_symbol'].notnull().any(): # gene_symbol already exists
                info['ensembl_id'] = info['gene_name']
            else:
       
                sample = str(info['gene_name'].dropna().iloc[0]).strip()
                
                if re.match('^ENS',sample):
                    info['ensembl_id'] = info['gene_name']
                    info['gene_symbol'] = None
                else:
                    info['gene_symbol'] = info['gene_name']
                    info['ensembl_id'] = None
            
        if 'gene_name' in info.columns:
            info.drop(columns=['gene_name'], inplace=True)
        
        
        
        expected_cols = ["log2fc", "logCPM", "pvalue", "padj", 'other_info', 'ensembl_id', 'gene_symbol']
        for col in expected_cols:
            if col not in info.columns:
                info[col] = None
        
        extra_columns = [col for col in info.columns if col not in expected_cols and col != 'experiment_id']
        
        if len(extra_columns) > 0:
            # Any extra input columns should be folded into `other_info`.
            extra_data = info[extra_columns].to_dict(orient='records')
            info.drop(columns=extra_columns, inplace=True)

            if 'other_info' in info.columns:
                # Standardize existing `other_info` values so dictionaries merge cleanly.
                existing_info = info['other_info'].apply(
                lambda x: json.loads(x) if isinstance(x, str) else (x if isinstance(x, dict) else {})
                )
        
                # Merge old and new dictionaries efficiently using a list comprehension
                info['other_info'] = [
                {**old, **new} if old else new 
                for old, new in zip(existing_info, extra_data)
                ]
            else:
            # If it didn't exist, safe to just assign it directly
                info['other_info'] = extra_data
            if 'other_info' not in info.columns:
                info['other_info'] = None

        
        
        info['experiment_id'] = id
        
        
        
        return info

    




    def initialize_gene_table(self, species):
        if species == 'human':
            url = 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'

        try:
            self.conn.execute('''
            INSERT INTO genes (symbol, id, ensembl_id, alias_symbol, prev_symbol, species)
            SELECT 
                csv_data.symbol, 
                CAST(regexp_replace(csv_data.hgnc_id, '^hgnc:', '', 'i') AS INTEGER) AS id, 
                csv_data.ensembl_gene_id, 
                string_split(UPPER(COALESCE(csv_data.alias_symbol, '')), '|') as alias_symbol, 
                string_split(UPPER(COALESCE(csv_data.prev_symbol, '')), '|') as prev_symbol, 
                'human' as species 
            FROM read_csv_auto(?) AS csv_data
            ''', (url,))
        except duckdb.ConstraintException as e:
            raise DuplicateGeneTableError(
                'Gene reference data appears to have already been loaded or the gene table has duplicate entries.'
            ) from e

    def _create_experiment(self, tool, date, file_path, data_signature, experiment_name=None, comparison_label=None):
        
        duplicate_ids=results = self.conn.execute(
            "SELECT experiment_id FROM experimental_data WHERE data_signature = ?", 
            (data_signature,)
            ).fetchall()
        if len(duplicate_ids) > 0:
            raise DuplicateExperimentError(
                f'Data is identical to data in experiment id: {duplicate_ids[0][0]}. Duplicate experiments are not allowed.'
            )


        
        result = self.conn.execute(
            "INSERT INTO experimental_data (tool, date, file, experiment_name, comparison_label, data_signature) VALUES (?, ?, ?, ?, ?, ?) RETURNING experiment_id",
            (tool, date, file_path, experiment_name, comparison_label, data_signature)
        ).fetchall()
        self.conn.commit()
        
        return result[0][0]

    def _insert_gene_results(self, insert_df):
        
        self.conn.execute('''
        INSERT INTO gene_results 
        (experiment_id, gene_symbol, ensembl_id, log2fc, logCPM, pvalue, padj, other_info) 
    
        SELECT DISTINCT ON (df.experiment_id, COALESCE(df.gene_symbol, df.ensembl_id))
            df.experiment_id, 
        
            -- 1. Grab the canonical symbol from the symbol match, or the ensembl match, or fallback to the input symbol
            COALESCE(ref_sym.symbol, ref_ens.symbol, df.gene_symbol) as gene_symbol, 
        
            -- 2. Grab the canonical ensembl_id from the ensembl match, or the symbol match, or fallback to the input ensembl_id
            COALESCE(ref_ens.ensembl_id, ref_sym.ensembl_id, df.ensembl_id) as ensembl_id,
        
            df.log2fc, 
            df.logCPM, 
            df.pvalue, 
            df.padj, 
            df.other_info
        
        FROM insert_df df
    
        -- JOIN 1: Try to resolve the user's gene_symbol column against official symbols, aliases, and past symbols
        LEFT JOIN genes ref_sym ON (
        UPPER(TRIM(df.gene_symbol)) = UPPER(ref_sym.symbol) OR 
        list_contains(ref_sym.alias_symbol, UPPER(TRIM(df.gene_symbol))) OR
        list_contains(ref_sym.prev_symbol, UPPER(TRIM(df.gene_symbol)))
        )
    
        -- JOIN 2: Try to resolve the user's ensembl_id column (stripping versions like .1) against official Ensembl IDs
        LEFT JOIN genes ref_ens ON (
        UPPER(TRIM(split_part(df.ensembl_id, '.', 1))) = UPPER(ref_ens.ensembl_id)
        )
    
        ORDER BY 
            df.experiment_id, 
            COALESCE(df.gene_symbol, df.ensembl_id), 
            ref_sym.symbol NULLS LAST, 
            ref_ens.ensembl_id NULLS LAST
        ''')

        self.conn.commit()
        
    def insert_to_database(self, info, tool=None, date=None, file_path=None, experiment_name=None, comparison_label=None):

        if experiment_name is None or comparison_label is None:
            raise ProcessingError('You must include an experiment name and comparison label for your experiment!')
        date, file_path = self._find_date_and_file(info)
        if date:
            date = pd.to_datetime(date).strftime('%Y-%m-%d')
        data_signature=self._get_data_signature(info)
        id = self._create_experiment(tool, date, file_path, data_signature, experiment_name, comparison_label)
        df = self._preprocess_df(info, id)
        self._insert_gene_results(df)
        

    def close(self):
        if self.conn:
            self.conn.close()
            self.conn = None
            
    def query(self, table, save_path=None, **filters):
        """Query a table using keyword filters.

        Filters are provided as kwargs like `pvalue__gt=0.05` or `gene_symbol='TP53'`.
        The method resolves filter column names to actual table columns; if a
        close (typo) match exists a ValueError is raised suggesting the correct
        name. If no close match exists and `other_info` is present, the filter
        will be applied against `other_info` JSON.
        Returns a pandas DataFrame with results.
        """
        table_names = [t[0] for t in self.tables]
        if table not in table_names:
           raise ProcessingError(f"Table '{table}' does not exist in the database. Available tables: {', '.join(table_names)}")

        if table == 'gene_results':
            actual_columns = [col[3] for col in self.gene_columns]
        elif table == 'experimental_data':
            actual_columns = [col[3] for col in self.experiment_columns]
        else:
            cols = self.conn.execute("SELECT * FROM information_schema.columns WHERE table_name=?", (table,)).fetchall()
            actual_columns = [col[3] for col in cols]

        clauses = []
        params = []
        for filter_key, filter_value in filters.items():
            filter_column, filter_operator = self._split_filter_key(filter_key)
            clause, clause_params = self._build_filter_clause(filter_column, filter_operator, filter_value, actual_columns, table)
            clauses.append(clause)
            params.append(clause_params)
        

        # Build a simple WHERE clause incrementally from normalized filter expressions.
        query = f'SELECT * FROM {table} WHERE 1=1'
        for clause in clauses:
            query += f' AND {clause}'

        
        df = self.conn.execute(query, params).df()
        
        if 'id' in df.columns:
            df = df.set_index('id')
        elif not save_path and 'experiment_id' in df.columns:
            df = df.set_index('experiment_id')

        if save_path:
            folder = os.path.dirname(save_path)
            if folder and not os.path.exists(folder):
                os.makedirs(folder)

            if save_path.endswith('.csv'):
                df.to_csv(save_path, index=False)
            elif save_path.endswith('.xlsx'):
                try:
                    import openpyxl
                    df.to_excel(save_path, index=False)
                except ImportError:
                    print("openpyxl library is required to save as Excel. Please install it using 'pip install openpyxl'.")
            elif save_path.endswith('.json'):
                df.to_json(save_path, orient='records')
            else:
                df.to_csv(save_path, index=False)

            print(f"Successfully saved results to {save_path}")

        print(f'found {len(df)} results matching the query.')
        return df
       
   
    
    def delete_gene_results(self, **filters):
        """Delete rows from `gene_results` matching provided filters.

        Filters use the same syntax as `query()`. This method resolves filters
        into SQL WHERE clauses and executes a DELETE. It will raise if no
        filters are provided to avoid accidental full-table deletes.
        """
        actual_columns = [col[3] for col in self.gene_columns]
        if not filters:
            raise ProcessingError("No filters provided. Don't delete everything by accident!")

        clauses = []
        params = []
        for filter_key, filter_value in filters.items():
            filter_column, filter_operator = self._split_filter_key(filter_key)
            clause, clause_params = self._build_filter_clause(filter_column, filter_operator, filter_value, actual_columns, 'gene_results')
            clauses.append(clause)
            params.append(clause_params)

        # Perform a filtered delete; this always requires filters to avoid full-table removal.
        query = 'DELETE FROM gene_results WHERE 1=1'
        for clause in clauses:
            query += f' AND {clause}'

        self.conn.execute(query, params)
        self.conn.commit()
        print('successfully deleted gene results')
    

        
        
        
       

    def delete_experiments(self, **filters):
        """Delete experiments (and associated gene_results) matching filters.

        Resolves filters into a WHERE clause against `experimental_data` and
        deletes matching experiments and any gene_results referencing them.
        """
        actual_columns = [col[3] for col in self.experiment_columns]
        if len(filters) == 0:
            raise ProcessingError("No filters provided for deletion.")

        clauses = []
        params = []
        for filter_key, filter_value in filters.items():
            filter_column, filter_operator = self._split_filter_key(filter_key)
            clause, clause_params = self._build_filter_clause(filter_column, filter_operator, filter_value, actual_columns, 'experimental_data')
            clauses.append(clause)
            params.append(clause_params)

        where_clause = " WHERE 1=1"
        for clause in clauses:
            where_clause += f' AND {clause}'

        try:
            self.conn.execute("BEGIN TRANSACTION")
            self.conn.execute(
                f"DELETE FROM gene_results WHERE experiment_id IN (SELECT experiment_id FROM experimental_data{where_clause})",
                params
            )
            self.conn.commit()
            self.conn.execute(f"DELETE FROM experimental_data{where_clause}", params)
            self.conn.commit()
            print('deleted experiments and associated gene results')
        except Exception as e:
            self.conn.execute("ROLLBACK")
            raise
       

        
    
    def execute_raw(self, query):
        self.conn.execute(query)

    

        