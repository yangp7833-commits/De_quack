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
2. Use `ingest()` to normalize and insert a DataFrame (preprocess_df).
3. `query()` and delete methods accept keyword filters, resolve them to
     concrete SQL clauses, and execute against DuckDB.
"""
import duckdb
from duckdb import SQLExpression, CaseExpression, ColumnExpression, ConstantExpression, FunctionExpression
import re
import os
import json
import uuid
from pathlib import Path
import hashlib
import sys
import difflib
import datetime
import pandas as pd
from .exceptions import ProcessingError, DuplicateExperimentError, DuplicateGeneTableError
from .utilities import gene_columns, ExperimentMetadata

class de_quackling:
    def __init__(self, db_path='SQL.duckdb'):
        """Initialize a de_quackling manager for a DuckDB-backed differential expression store.

        Args:
            db_path: path to the DuckDB file to use or create.
        """
        self.db_path = db_path
        self.conn = None
        
        
        
    def __enter__(self):
        """Open the DuckDB connection and initialize required tables.

        This method is used by the context manager protocol and creates the
        core schema (`experimental_data`, `gene_results`, `genes`) plus indexes
        and sequences if they do not already exist.
        """
        self.conn = duckdb.connect(self.db_path)

        # Maintain a sequence for experiment primary keys.
        self.conn.execute('''CREATE SEQUENCE IF NOT EXISTS seq_experiment_id START 1''')

        # Create tables if they don't already exist.
        self.conn.execute('''CREATE TABLE IF NOT EXISTS experimental_data
                            (experiment_id INTEGER PRIMARY KEY DEFAULT nextval('seq_experiment_id'), 
                             model VARCHAR, date VARCHAR, file VARCHAR, 
                             experiment_name VARCHAR, contrast VARCHAR, annotation_version VARCHAR, normalization VARCHAR, other_info VARCHAR, data_signature VARCHAR)''')
        
        
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
        
        self.conn.execute('''CREATE TABLE IF NOT EXISTS gene_columns 
                        (column_name VARCHAR, alias_names VARCHAR[]
                        )''')
        self.conn.execute('CREATE UNIQUE INDEX IF NOT EXISTS idx_gene_columns_name ON gene_columns(column_name)')
        self.conn.execute('''INSERT OR REPLACE INTO gene_columns (column_name, alias_names) 
                        VALUES ('gene_symbol', ?), ('ensembl_id', ?), ('pvalue', ?), ('padj', ?), ('log2fc', ?), ('base_mean', ?), ('logCPM', ?), ('gene_name', ?)''',
                        (gene_columns['gene_symbol'], gene_columns['ensembl_id'], gene_columns['pvalue'], gene_columns['padj'], gene_columns['log2fc'], gene_columns['base_mean'], gene_columns['logCPM'], gene_columns['gene_name']))

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
        """Close the DuckDB connection when leaving the context manager."""
        if self.conn:
            self.conn.close()
    
    def connect(self):
        """Ensure the database connection is open and return the manager instance."""
        if not self.conn:
            self.__enter__()
        return self



    def _get_data_signature(self):
        sample_data=self.conn.execute('SELECT * FROM (SELECT * FROM preprocessed_data LIMIT 100)').fetchall()
    
        # Generate a unique string (hash) from that data
        return hashlib.md5(str(sample_data).encode()).hexdigest()

   
    
        

    def _create_temp_view(self): 
        from duckdb import SQLExpression

    
        columns_info = self.conn.execute(f'''
            SELECT 
                t.column_name AS incoming_col, 
                g.column_name AS real_col
            FROM duckdb_columns() AS t
            LEFT JOIN gene_columns AS g
                ON list_contains(g.alias_names, lower(t.column_name))
            WHERE t.table_name = 'preprocessed_data'
            ORDER BY t.column_index
        ''').fetchall()

        if not columns_info:
            raise ProcessingError("The registered table contains no columns.")
      
        temp_view = self.conn.table('preprocessed_data')

        mapped_incoming = []
        all_columns = []
        columns = []
        expressions=[]
        numeric_columns = {'log2fc', 'logCPM', 'pvalue', 'padj', 'base_mean'}

        for incoming_col, real_col in columns_info:
            all_columns.append(incoming_col)
            # real_col may be None if no alias matched
            if real_col == 'gene_name':
                # sample a real cell value from this incoming column (avoid placeholder misuse)
                sample = temp_view.select(ColumnExpression(incoming_col)).limit(1).fetchall()
                sample_val = sample[0] if sample is not None else None
                if sample_val is not None and re.match(r'^ENS', str(sample_val)):
                    mapped_incoming.append(incoming_col)
                    expressions.append(ColumnExpression(incoming_col).alias('ensembl_id'))
                    columns.append('ensembl_id')
                elif 'gene_symbol' not in [col[0] for col in columns_info]:
                    mapped_incoming.append(incoming_col)
                    expressions.append(ColumnExpression(incoming_col).alias('gene_symbol'))
                    columns.append('gene_symbol')
                else:
                    pass
            elif real_col is not None:
                if real_col != 'base_mean':
                    mapped_incoming.append(incoming_col)
                expressions.append(ColumnExpression(incoming_col).alias(real_col))
                columns.append(real_col)

        if not mapped_incoming:
            raise ProcessingError(
                'No recognizable gene result columns were found in the input. '
                'Verify the header names against gene_columns aliases.'
            )
    
        extra_columns = [col for col in all_columns if col not in mapped_incoming]
        if extra_columns:
            struct_items=[]
            for column in extra_columns:
                struct_items.append(ConstantExpression(column))
                struct_items.append(ColumnExpression(column))
            json_expr = FunctionExpression('json_object', *struct_items)
            expressions.append(json_expr.alias('other_info'))
        else:
            expressions.append(ConstantExpression(None).alias('other_info'))
        
        
        
        if 'gene_symbol' not in columns:
            expressions.append(SQLExpression('NULL').alias('gene_symbol'))
        if 'ensembl_id' not in columns:
            expressions.append(SQLExpression('NULL').alias('ensembl_id'))
        if 'logCPM' not in columns:
            expressions.append(SQLExpression('NULL').alias('logCPM'))

        if not mapped_incoming:
            raise ProcessingError(
                'No recognizable gene result columns were found in the input. '
                'Verify the header names against gene_columns aliases.'
            )
        temp_view=temp_view.select(*expressions)
        return temp_view




    def _preprocess(self, info):

        if isinstance(info, str) and os.path.isfile(os.path.abspath(info)):
            info=os.path.abspath(info)
            if info.lower().endswith('.parquet'):
                self.conn.execute('DROP VIEW IF EXISTS preprocessed_data')
                self.conn.read_parquet(info, header=True).create_view('preprocessed_data')
                return
            else:
                self.conn.execute('DROP VIEW IF EXISTS preprocessed_data')
                self.conn.from_csv_auto(info, header=True).create_view("preprocessed_data")
                return
        if isinstance(info, str) and not os.path.isfile(os.path.abspath(info)):
            raise ProcessingError(f"{os.path.abspath(info)} is not a valid file path.")
        

        class_name = info.__class__.__name__.lower()
        self.conn.execute('DROP VIEW IF EXISTS preprocessed_data')

        if isinstance(info, de_arrow):
            self.conn.execute('CREATE TEMP VIEW preprocessed_data AS SELECT * FROM info')
            return


        if class_name == 'dataframe':
            self.conn.execute('CREATE TEMP VIEW preprocessed_data AS SELECT * FROM info')
            return

        if class_name == 'Table':
            self.conn.execute('CREATE TEMP VIEW preprocessed_data AS SELECT * FROM info')
            return
        

        raise ProcessingError(f'type {type(info)} is not supported. Provide a pandas DataFrame, list of dicts, or path to a CSV/TSV/Parquet file.')


    

    def initialize_gene_table(self, species):
        """Load species-specific reference gene annotation data into the database.

        Currently supports `human` by downloading HGNC gene metadata and
        populating the `genes` reference table for symbol and Ensembl matching.
        """
        if species == 'human':
            url = 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'
        else:
            raise ProcessingError(f"Unsupported species '{species}'. Only 'human' is currently supported.")

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
            self.conn.commit()

        except duckdb.ConstraintException as e:
            raise DuplicateGeneTableError(
                'Gene reference data appears to have already been loaded or the gene table has duplicate entries.'
            ) from e

    def _create_experiment(self,metadata_fields, data_signature, other_info):
        duplicate_ids = self.conn.execute(
            "SELECT experiment_id FROM experimental_data WHERE data_signature = ?", 
            (data_signature,)
        ).fetchall()
        #if len(duplicate_ids) > 0:
            #raise DuplicateExperimentError(
                #f'Data is identical to data in experiment id: {duplicate_ids[0][0]}. Duplicate experiments are not allowed.'
            #)
        result = self.conn.execute(
            "INSERT INTO experimental_data (model, date, file, experiment_name, contrast, annotation_version, normalization, other_info, data_signature) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?) RETURNING experiment_id",
            (metadata_fields.get('model'), metadata_fields.get('date'), metadata_fields.get('file'), metadata_fields.get('experiment_name'), metadata_fields.get('contrast'), metadata_fields.get('annotation_version'), metadata_fields.get('normalization'), other_info, data_signature)
        ).fetchall()
        self.conn.commit()
        
        return result[0][0]

    def _insert_gene_results(self, experiment_id, view):


        self.conn.execute('''
        INSERT INTO gene_results (experiment_id, gene_symbol, ensembl_id, log2fc, logCPM, pvalue, padj, other_info)
        SELECT ? AS experiment_id,
            COALESCE(ens.symbol, sym.symbol, e.gene_symbol) AS gene_symbol,
            COALESCE(ens.ensembl_id, sym.ensembl_id, e.ensembl_id) AS ensembl_id,
            COALESCE(TRY_CAST(e.log2fc AS DOUBLE), NULL),
            COALESCE(TRY_CAST(e.logCPM AS DOUBLE), log2(e.base_mean + 1), NULL),
            COALESCE(TRY_CAST(e.pvalue AS DOUBLE), NULL),
            COALESCE(TRY_CAST(e.padj AS DOUBLE), NULL),
            e.other_info
        FROM view as e
        LEFT JOIN genes ens ON  (
            e.ensembl_id IS NOT NULL AND ens.ensembl_id = UPPER(TRIM(e.ensembl_id))
        )

        LEFT JOIN genes sym ON (
            e.gene_symbol IS NOT NULL AND sym.symbol = UPPER(TRIM(e.gene_symbol)) OR
            e.gene_symbol IS NOT NULL AND sym.prev_symbol IS NOT NULL AND list_contains(sym.prev_symbol, e.gene_symbol)
            )
        
        
        ''', (experiment_id,))

        self.conn.commit()
        
    def ingest(self, info, metadata: dict, **config_columns):
        """Normalize and ingest differential expression results into DuckDB.

        Args:
            metadata: A dictionary containing the experiment metadata.
            **config_columns: additional configuration columns for the ingestion process.
        

            experiment_name: human-readable experiment name.
            date: optional experiment date; auto-detected from file path or current time.
            tool: optional source tool name.
            **config_columns: additional configuration columns for the ingestion process.

        The method validates required metadata, calculates a data signature,
        prevents duplicate experiments, and loads normalized gene results into
        `gene_results` with canonical symbol/Ensembl mapping.
        """

        metadata_fields, other_info = ExperimentMetadata().to_dict(metadata, info)

        self._preprocess(info)
        data_signature = self._get_data_signature()
        id = self._create_experiment(metadata_fields, data_signature, other_info)
        view=self._create_temp_view()
        self._insert_gene_results(id, view)
    

    def close(self):
        """Close the active DuckDB connection and reset internal state."""
        if self.conn:
            self.conn.close()
            self.conn = None
    



            
       

        
    
    def execute_raw(self, query, values=None):
        return self.conn.execute(query, (values,))



    

        