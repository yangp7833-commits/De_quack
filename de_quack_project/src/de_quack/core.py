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
import nanoarrow as na
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
from .exceptions import ProcessingError, DuplicateExperimentError, DuplicateGeneTableError, DeQuackError
from .utilities import gene_columns, ExperimentMetadata, CORE_QUERIES, gene_mapping
import time
_gene_mapping_queries = gene_mapping
core_queries = CORE_QUERIES
experiment_columns=['experiment_id', 'model', 'date', 'file', 'experiment_name', 'contrast', 'annotation_version', 'normalization', 'extra_info']

_GENE_ALIAS_TO_COLUMN = {
    alias.lower(): canonical
    for canonical, aliases in gene_columns.items()
    for alias in aliases + [canonical]
}

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
        tables = self.conn.execute("SELECT table_name FROM information_schema.tables WHERE table_schema='main'").fetchall()
        self.tables = [table[0] for table in tables]
        if 'gene_results' not in self.tables:
            self._initialize_tables()

        self.gene_columns = self.conn.execute("SELECT * FROM information_schema.columns WHERE table_name='gene_results'").fetchall()
        self.experiment_columns = self.conn.execute("SELECT * FROM information_schema.columns WHERE table_name='experimental_data'").fetchall()
        return self
    
    def _initialize_tables(self):
        self.conn.execute('''CREATE SEQUENCE IF NOT EXISTS seq_experiment_id START 1''')

        # Create tables if they don't already exist.
        self.conn.execute('''CREATE TABLE IF NOT EXISTS experimental_data
                            (experiment_id INTEGER PRIMARY KEY DEFAULT nextval('seq_experiment_id'), 
                             model VARCHAR, date VARCHAR, file VARCHAR, 
                             experiment_name VARCHAR, contrast VARCHAR, annotation_version VARCHAR, normalization VARCHAR, other_info VARCHAR, data_signature VARCHAR)''')
        
        
        
        self.conn.execute('''CREATE TABLE IF NOT EXISTS gene_results
                            (experiment_id INTEGER,
                            gene_symbol VARCHAR,
                            ensembl_id VARCHAR,
                            log2fc DOUBLE,
                            logCPM DOUBLE,
                            pvalue DOUBLE,
                            padj DOUBLE,
                            stat DOUBLE,
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
    
    def __exit__(self, exc_type, exc_value, traceback):
        """Close the DuckDB connection when leaving the context manager."""
        if self.conn:
            self.conn.close()

    def fast_connect(self):
        self.conn = duckdb.connect(self.db_path)
        return self
    
    def connect(self):
        """Ensure the database connection is open and return the manager instance."""
        if not self.conn:
            self.__enter__()
        return self



    def _get_data_signature(self):
        sample_data=self.conn.execute('SELECT * FROM (SELECT * FROM preprocessed_data LIMIT 100)').fetchall()
    
        # Generate a unique string (hash) from that data
        return hashlib.md5(str(sample_data).encode()).hexdigest()

   
    
        

    def _create_temp_view(self, columns = {}): 
        from duckdb import SQLExpression

        temp_view = self.conn.table('preprocessed_data')
        columns_info = [(column_name, _GENE_ALIAS_TO_COLUMN.get(column_name.lower())) for column_name in temp_view.columns]
        for key, value in columns.items():
            if key in temp_view.columns and value in gene_columns:
                if (key, value) not in columns_info:
                    columns_info.append((key, value))
            else:
                raise ProcessingError(f"Invalid column mapping: {key} -> {value}. Ensure the column exists in the data and the value is a valid gene column alias.")


        if not columns_info:
            raise ProcessingError("The registered table contains no columns.")

        has_gene_symbol = any(real_col == 'gene_symbol' for _, real_col in columns_info)

        mapped_incoming = []
        all_columns = []
        columns = []
        expressions=[]

        for incoming_col, real_col in columns_info:
            all_columns.append(incoming_col)
            if real_col == 'gene_name':
                sample = temp_view.select(ColumnExpression(incoming_col)).limit(1).fetchone()
                sample_val = sample[0] if sample is not None else None
                if sample_val is not None and re.match(r'^ENS', str(sample_val)):
                    mapped_incoming.append(incoming_col)
                    expressions.append(ColumnExpression(incoming_col).alias('ensembl_id'))
                    columns.append('ensembl_id')
                elif not has_gene_symbol:
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
        for col in ['gene_symbol', 'ensembl_id', 'log2fc', 'padj', 'pvalue', 'stat', 'logCPM']:
            if col not in columns:
                expressions.append(ConstantExpression(None).alias(col))
                columns.append(col)
    
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
        if 'base_mean' not in columns:
            expressions.append(SQLExpression('NULL').alias('base_mean'))

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

        if class_name in {'de_arrow', 'de_arrows'}:
            self.conn.register('info', info._table)
            self.conn.execute('CREATE TEMP VIEW preprocessed_data AS SELECT * FROM info')
            return

        try:
            self.conn.register('info', info)
        except Exception:
            pass

        if class_name == 'dataframe':
            self.conn.execute('CREATE TEMP VIEW preprocessed_data AS SELECT * FROM info')
            return

        if class_name == 'table':
            self.conn.execute('CREATE TEMP VIEW preprocessed_data AS SELECT * FROM info')
            return
        

        raise ProcessingError(f'type {type(info)} is not supported. Provide a pandas DataFrame, list of dicts, or path to a CSV/TSV/Parquet file.')
    
    def initialize_gene_table(self, species = 'human'):
        """Initialize the reference gene table for a given species.

        Currently supports `human` by downloading HGNC gene metadata and
        populating the `genes` reference table for symbol and Ensembl matching.
        """
        if species.lower() == 'human':
            self._insert_human_genes()
        elif species.lower() == 'mouse':
            self._insert_mouse_genes()
        else:
            raise ProcessingError(f"Species '{species}' is not supported. Only 'human' is currently supported.")
    
    def _insert_mouse_genes(self):
        """Load species-specific reference gene annotation data into the database.

        Currently supports `mouse` by downloading MGI gene metadata and
        populating the `genes` reference table for symbol and Ensembl matching."""

       
        try:
            self.conn.execute(_gene_mapping_queries['mouse_genes'])
            self.conn.commit()
        except duckdb.ConstraintException as e:
            raise DuplicateGeneTableError(
                'Gene reference data appears to have already been loaded or the gene table has duplicate entries.'
            ) from e

    

    def _insert_human_genes(self):
        """Load species-specific reference gene annotation data into the database.

        Currently supports `human` by downloading HGNC gene metadata and
        populating the `genes` reference table for symbol and Ensembl matching."""

        try:
            self.conn.execute(_gene_mapping_queries['human_genes'])
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

    def _insert_gene_results(self, experiment_id, view, species = None):


        self.conn.execute('''
        INSERT INTO gene_results (experiment_id, gene_symbol, ensembl_id, log2fc, logCPM, pvalue, padj, stat, other_info)
        SELECT $1 AS experiment_id,
            COALESCE(ens.symbol, sym.symbol, e.gene_symbol) AS gene_symbol,
            COALESCE(ens.ensembl_id, sym.ensembl_id, e.ensembl_id) AS ensembl_id,
            COALESCE(TRY_CAST(e.log2fc AS DOUBLE), NULL) AS log2fc,
            COALESCE(TRY_CAST(e.logCPM AS DOUBLE), log2(e.base_mean + 1), NULL) AS logCPM,
            COALESCE(TRY_CAST(e.pvalue AS DOUBLE), NULL) AS pvalue,
            COALESCE(TRY_CAST(e.padj AS DOUBLE), NULL) AS padj,
            COALESCE(TRY_CAST(e.stat AS DOUBLE), NULL) AS stat,
            e.other_info
        FROM view as e
        LEFT JOIN genes ens ON  (
            e.ensembl_id IS NOT NULL AND ens.ensembl_id = UPPER(TRIM(e.ensembl_id))
        )

        LEFT JOIN genes sym ON 
            (e.gene_symbol IS NOT NULL AND sym.symbol = e.gene_symbol) OR
            e.gene_symbol IS NOT NULL AND sym.prev_symbol IS NOT NULL AND list_contains(sym.prev_symbol, e.gene_symbol)
            
        
        WHERE (sym.species IS NULL OR sym.species = $2) AND (ens.species IS NULL OR ens.species = $2)
        
        
        ''', (experiment_id, species))

        self.conn.commit()
        
    def ingest(self, info, metadata: dict = {}, species = 'human', columns = {}, **kwargs):
        """Normalize and ingest differential expression results into DuckDB.

        Args:
            metadata: A dictionary containing the experiment metadata.
            config_columns: A dictionary of additional configuration columns for the ingestion process.
            **kwargs: Additional keyword arguments.

            experiment_name: human-readable experiment name.
            date: optional experiment date; auto-detected from file path or current time.
            tool: optional source tool name.
            **config_columns: additional configuration columns for the ingestion process.

        The method validates required metadata, calculates a data signature,
        prevents duplicate experiments, and loads normalized gene results into
        `gene_results` with canonical symbol/Ensembl mapping.
        """
        if kwargs:
            metadata.update(kwargs)
        metadata_fields, other_info = ExperimentMetadata().to_dict(metadata, info)
        self._preprocess(info) 

        data_signature = self._get_data_signature()

        id = self._create_experiment(metadata_fields, data_signature, other_info)

        view=self._create_temp_view(columns)

        self._insert_gene_results(id, view, species)
        

    def get_experiment(self, id = None, name = None, model = None, annotation_version = None,  normalization = None, date = None, contrast = None, file = None):
        try:
            date = datetime.datetime.strptime(date, '%Y-%m-%d').date() if date else None
        except ValueError:
            raise ProcessingError(f"Invalid date format: {date}. Expected format is YYYY-MM-DD.")
        rel = self.conn.execute(core_queries['get_experiment'], [id, date, model, file, name, contrast, annotation_version, normalization]).fetchall()
        metadata_list = [dict(zip(experiment_columns, row)) for row in rel]
        metadata_fields, ids = _to_metadata(metadata_list)
        self.conn.execute(f'CREATE OR REPLACE TEMP TABLE ids AS SELECT UNNEST({ids}) AS id')
        df = self.conn.execute('SELECT * exclude(ids.id) FROM gene_results g JOIN ids ON g.experiment_id = ids.id').pl()
        return _get_de_arrows(df, metadata_fields, ids)
        
    
    def get_significant_genes(self, log2fc = 1, padj = 0.05, logCPM = 1, id = None):
        df = self.conn.execute('SELECT * FROM gene_results WHERE abs(log2fc) > $1 AND padj < $2 AND logCPM > $3 AND (experiment_id = $4 OR $4 IS NULL)', [log2fc, padj, logCPM, id]).pl()
        return self._polars_to_de_arrows(df)

    def get_upregulated(self, log2fc = 1, padj = 0.05, logCPM = 1, id = None):
        df = self.conn.execute('SELECT * FROM gene_results WHERE log2fc > $1 AND padj < $2 AND logCPM > $3 AND (experiment_id = $4 OR $4 IS NULL)', [log2fc, padj, logCPM, id]).pl()
        arrow = self._polars_to_de_arrows(df)
        return arrow
        
    
    def get_downregulated(self, log2fc = -1, padj = 0.05, logCPM = 1, id = None):
        df = self.conn.execute('SELECT * FROM gene_results WHERE log2fc < $1 AND padj < $2 AND logCPM > $3 AND (experiment_id = $4 OR $4 IS NULL)', [log2fc, padj, logCPM, id]).pl()
        return self._polars_to_de_arrows(df)

    def get_gene(self, gene_symbol = None, ensembl_id = None, id = None):
        rel = self.conn.execute('SELECT * FROM gene_results WHERE (gene_symbol = $1 OR $1 IS NULL) AND (ensembl_id = $2 OR $2 IS NULL) AND (experiment_id = $3 OR $3 IS NULL)', [gene_symbol, ensembl_id, id]).pl()
        return self._polars_to_de_arrows(rel)

    def delete_experiment(self, id = None, name = None, model = None, annotation_version = None,  normalization = None, date = None, contrast = None, file = None):
        try:
            date = datetime.datetime.strptime(date, '%Y-%m-%d').date() if date else None
        except ValueError:
            raise ProcessingError(f"Invalid date format: {date}. Expected format is YYYY-MM-DD.")
        try:
            self.conn.begin()
            ids = self.conn.execute(core_queries['find_delete_experiment'], (id, date, model, file, name, contrast, annotation_version, normalization)).fetchall()
            ids = [row[0] for row in ids]
            self.conn.execute(core_queries['delete_gene_results'], (ids,))
            self.conn.commit()
            self.conn.execute(core_queries['delete_experiment'], (ids,))
        finally:
            self.conn.commit()
    
    def _polars_to_de_arrows(self, df):
        if len(df.columns) == 0:
            return self._get_de_arrows(df, {}, [])
        ids = df['experiment_id'].unique().to_list()
        meta_rel = self.conn.execute(
            '''
            SELECT
                experiment_id,
                model,
                date,
                file,
                experiment_name,
                contrast,
                annotation_version,
                normalization,
                other_info AS extra_info
            FROM experimental_data
            WHERE experiment_id = ANY(?)
            ''',
            [ids]
        )
        metadata_list = [dict(zip(experiment_columns, row)) for row in meta_rel.fetchall()]
        metadata_fields, ids = _to_metadata(metadata_list)
        return _get_de_arrows(df, metadata_fields, ids)
        
            



    
        
        

    

    def close(self):
        """Close the active DuckDB connection and reset internal state."""
        if self.conn:
            self.conn.close()
            self.conn = None



    def execute_raw(self, query):
        return self.conn.execute(query)
    
def _get_de_arrows(arrows, metadata, ids):
    from .arrow import DeArrow, DeArrows
    if len(ids) < 2:
        return DeArrow._from_arrow(arrows, metadata, ids)
    return DeArrows._from_arrow(arrows, metadata, ids)

def _to_metadata(metadata_list):
    for meta in metadata_list:
        if meta.get('extra_info'):
            try:
                extra = json.loads(meta['extra_info'])
                for key, value in extra.items():
                    meta[key] = value
            except json.JSONDecodeError:
                pass
        meta.pop('extra_info', None)
    ids = [meta.pop('experiment_id') for meta in metadata_list]
    metadata_fields = {exp_id: meta for meta, exp_id in zip(metadata_list, ids)}
    return metadata_fields, ids





    

        