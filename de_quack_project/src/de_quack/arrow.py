import nanoarrow as na
from .core import de_quackling
from .utilities import ExperimentMetadata, DE_ARROW_QUERIES, DE_ARROWS_QUERIES, gene_columns, gene_mapping
import duckdb
from .exceptions import DeQuackError, ProcessingError, DuplicateExperimentError, DuplicateGeneTableError
import json
import uuid
import re
import os
import sys
from importlib import resources
import hashlib
import time
import polars as pl
_de_arrow_queries=DE_ARROW_QUERIES
_de_arrows_queries=DE_ARROWS_QUERIES
_gene_mapping_queries = gene_mapping
_GENE_ALIAS_TO_COLUMN = {
    alias.lower(): canonical
    for canonical, aliases in gene_columns.items()
    for alias in aliases + [canonical]
}
_gene_columns = ['gene_symbol', 'ensembl_id', 'log2fc', 'logCPM', 'pvalue', 'padj', 'stat']
class de_arrow:
    
    def __init__(self, info, id = None, metadata=None,  heal_genes = False, species = None, **fields):
        
        if metadata is None:
            metadata = fields
        if len(metadata) == 0:
            raise ProcessingError('No metadata provided in de_arrow object')
        
        if id is None:
            id = 1
        
        self._table, self.experiment_metadata = self.__class__._to_de_arrow(info, metadata, id, heal_genes = heal_genes, species = None)
        self.name = _get_experiment_attribute(self.experiment_metadata, 'experiment_name', id)
        self.annotation_version = _get_experiment_attribute(self.experiment_metadata, 'annotation_version', id)
        self.contrast = _get_experiment_attribute(self.experiment_metadata, 'contrast', id)
        self.file = _get_experiment_attribute(self.experiment_metadata, 'file', id)
        self.date = _get_experiment_attribute(self.experiment_metadata, 'date', id)
        self.model = _get_experiment_attribute(self.experiment_metadata, 'model', id)
        self.id = id

        
    
    def __getattr__(self, name):
        # Delegate attribute access to the internal nanoarrow table without triggering recursion
        table = object.__getattribute__(self, '_table')
        return getattr(table, name)
    
    def __setattr__(self, name, value):
        if name in ('_table', 'experiment_metadata', 'name', 'annotation_version', 'contrast', 'file', 'date', 'id', 'model'):
            object.__setattr__(self, name, value)
        else:
            table = object.__getattribute__(self, '_table')
            setattr(table, name, value)
    
    def __getitem__(self, key):
        return self._table[key]
       

    def __len__(self):
        return len(self._table)
    
    def __repr__(self):
        return repr(self._table)
    
    def __str__(self):
        return(str(self._table))
    

    
    @classmethod
    def _to_de_arrow(self, info, metadata, experiment_id=None, heal_genes=False, species = None):
        df = _to_polars_table(info)
        metadata_fields, extra_info=ExperimentMetadata().to_dict(metadata, info)
        for key, value in json.loads(extra_info).items():
            metadata_fields[key]=value
        metadata_fields={experiment_id: metadata_fields}
        df = _order_columns(df)
        if heal_genes == True:
            if species is None:
                species = 'human'
            start_time = time.perf_counter()
            df = _heal_genes(df, species)
            end_time = time.perf_counter()
            elapsed_time = end_time - start_time
            print(f"Healing genes took {elapsed_time:.2f} seconds.")
        df = df.insert_column(0, pl.lit(experiment_id).alias('experiment_id'))
        return df, metadata_fields


   

    @classmethod
    def _from_arrow(cls, arrow_array, metadata, id, columns = None):
        """
        Convert a nanoarrow array back into a de_arrow object.
        
        Parameters
        ----------
        arrow_array : nanoarrow.Array
            A nanoarrow array to convert
        metadata : dict
            A nested dictionary containing metadata {ID: {field: value}}
        
        Returns
        -------
        de_arrow
            A de_arrow object with the provided array and metadata
        """
        instance = cls.__new__(cls)
        instance._table = arrow_array
        instance.experiment_metadata = metadata
        instance.name = _get_experiment_attribute(metadata, 'experiment_name', id)
        instance.annotation_version = _get_experiment_attribute(metadata, 'annotation_version', id)
        instance.contrast = _get_experiment_attribute(metadata, 'contrast', id)
        instance.file = _get_experiment_attribute(metadata, 'file', id)
        instance.date = _get_experiment_attribute(metadata, 'date', id)
        instance.model = _get_experiment_attribute(metadata, 'model', id)
        instance.id = id
        return instance
    
    
    def get_significant_genes(self, log2fc=1, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        de=de_quackling(get_unique_conn()).fast_connect()
        try:
            arrow_table=self._table
            rel=de.conn.sql(_de_arrow_queries['get_important_genes'], params={'log2fc':log2fc, 'pvalue':pvalue, 'logCPM': logCPM})
            arrow_table=na.Array(rel)
            arrow_table=self._from_arrow(arrow_table, self.experiment_metadata, self.id)
            return arrow_table
        finally:
            de.conn.close()
    
    def get_downregulated(self, log2fc=0, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        de=de_quackling(get_unique_conn()).fast_connect()
        try:
            arrow_table=self._table
            rel=de.conn.sql(_de_arrow_queries['get_downregulated'], params={'log2fc':log2fc, 'pvalue':pvalue, 'logCPM': logCPM})
            arrow_table=self._from_arrow(arrow_table, self.experiment_metadata, self.id)
            return arrow_table
        finally:
            de.conn.close()
    
    def get_upregulated(self, log2fc=0, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        de=de_quackling(get_unique_conn()).fast_connect()
        try:
            arrow_table=self._table
            rel=de.conn.sql(_de_arrow_queries['get_upregulated'], params={'log2fc':log2fc, 'pvalue':pvalue, 'logCPM': logCPM})
            arrow_table=na.Array(rel)
            arrow_table=self._from_arrow(arrow_table, self.experiment_metadata, self.id)
            return arrow_table
        finally:
            de.conn.close()
    
    def set_id(self, id):
        de=de_quackling(get_unique_conn()).fast_connect()
        try:
            arrow_table=self._table
            rel=de.conn.sql(_de_arrow_queries['set_arrow_id'], params={'id':id})
            arrow_table=na.Array(rel)
            old_id = self.id
            updated_metadata = {id: self.experiment_metadata.get(old_id, {})}
            arrow_table=self._from_arrow(arrow_table, updated_metadata, id)
            return arrow_table
        finally:
            de.conn.close()
    
    def add_experiment(self, data, metadata = None, id = None):
        if isinstance(data, de_arrow):
            return self._add_experiment_arrow(data, metadata, id)
        if isinstance(data, de_arrows):
            if metadata is not None:
                if not isinstance(metadata, list) or not isinstance(metadata[0], dict):
                    raise DeQuackError('Metadata must be a list of dictionaries for de_arrows objects!')
                if len(id) != len(data.id):
                    raise DeQuackError('Length of id list does not match length of incoming de_arrows object')
            return self._add_experiment_arrows(data, metadata, id)
        return self._add_experiment_data(data, metadata, id)

    
    def _add_experiment_arrows(self, data, metadata = None, ids = None):
        de = de_quackling(get_unique_conn()).fast_connect()
        try:
            arrow_table = self._table
            data = data.set_id(ids)
            if ids is None:
                ids = data.id
            ids.append(self.id)
            _check_ids(ids)
            self_metadata = self.experiment_metadata
            if metadata is None:
                new_metadata = data.experiment_metadata
                for id, value in zip(ids, new_metadata.values()):
                    self_metadata[id] = value
            else:
                new_metadata = metadata
                for meta, id in zip(ids, new_metadata):
                    m, extra_info = ExperimentMetadata.to_dict(new_metadata)
                    for key, value in extra_info.items():
                        m[key] = value
                    self_metadata[n] = m
            
            arrow_table = de.conn.sql("SELECT * FROM arrow_table UNION ALL SELECT * FROM data")
            arrow_table = na.Array(arrow_table)
            ids = [id for id in self_metadata.keys()]
            return de_arrows._from_arrow(arrow_table, self_metadata, ids)
        finally:
            de.conn.close()
        
        
    def _add_experiment_arrow(self, data, metadata = None, id = None):
        arrow_table = self._table
        if not isinstance(id, int) and id is not None:
            raise DeQuackError('id has to be an integer')
        de = de_quackling(get_unique_conn()).fast_connect()
        try:
            new_metadata = ExperimentMetadata.to_dict(metadata) if metadata is not None else data.experiment_metadata
            id = id if id is not None else data.id
            if id == self.id:
                raise DeQuackError('Duplicate ids found')
            metadata = self.experiment_metadata
            arrow_table = de.conn.sql("SELECT * FROM arrow_table UNION ALL SELECT * REPLACE ($id AS experiment_id) FROM data", params = {'id':id})
            arrow_table = na.Array(arrow_table)
            metadata[id] = list(new_metadata.values())[0]
            ids = [self.id, id]
            return de_arrows._from_arrow(arrow_table, metadata, ids)
        finally:
            de.conn.close()

    def _add_experiment_data(self, data, metadata = None, id = None):
        arrow_table = self._table
        if metadata is None:
            raise DeQuackError('No metadata provided for the data')
        if id is None:
            raise DeQuackError('No id provided')
        de = de_quackling(get_unique_conn()).connect()
        try:
            new_metadata, extra_info=ExperimentMetadata().to_dict(metadata, data)
            for key, value in json.loads(extra_info).items():
                new_metadata[key]=value
            new_metadata={id: new_metadata}
            de._preprocess(data)
            de_arrow_insertion_view = de._create_temp_view()
            data = de.conn.sql(_de_arrow_queries['insert_to_de_arrow'], params={'id':id})
            if id == self.id:
                raise DeQuackError('Duplicate ids found')
            metadata = self.experiment_metadata
            arrow_table = de.conn.sql("SELECT * FROM arrow_table UNION ALL SELECT * FROM data")
            arrow_table = na.Array(arrow_table)
            metadata[id] = new_metadata[id]
            ids = [self.id, id]
            return de_arrows._from_arrow(arrow_table, metadata, ids)
        finally:
            de.conn.close()
    
    def get_gene(self, gene_symbol = None, ensembl_id = None):
        required_columns = {'gene_symbol', 'ensembl_id'}
        _check_columns(required_columns, self.columns)
        de = de_quackling(get_unique_conn()).fast_connect()
        arrow_table = self._table
        rel = de.conn.sql(_de_arrow_queries['get_gene'], 
        params = {
        'gene_symbol': gene_symbol,
        'ensembl_id': ensembl_id,
        })
        arrow_table = na.Array(rel)
        de.conn.close()
        return self._from_arrow(arrow_table, self.experiment_metadata, self.id)

    def insert(self, file):
        required_columns = {'padj', 'pvalue', 'log2fc', 'gene_symbol', 'ensembl_id', 'logCPM', 'stat', 'other_info', 'experiment_id'}
        _check_columns(required_columns, self.columns)
        if not isinstance(file, str):
            raise TypeError(f'file must be a string, not {type(file)}')
        arrow_table = self._table
        metadata = self.experiment_metadata   
        m, other_info = ExperimentMetadata().to_dict(metadata[self.id], None)
        de = de_quackling(file).connect()
        try:
            sample_data = de.conn.execute('SELECT * FROM (SELECT * FROM arrow_table LIMIT 100)').fetchall()
            data_sig = hashlib.md5(str(sample_data).encode()).hexdigest()
            id = de._create_experiment(m, data_sig, other_info)
            de.conn.sql('INSERT INTO gene_results (experiment_id, gene_symbol, ensembl_id, log2fc, logCPM, pvalue, padj, stat, other_info) SELECT * REPLACE($id AS experiment_id) FROM arrow_table', params = {'id': id})
            de.conn.commit()
            return self
        finally:
            de.conn.close()

    
    def df(self):
        try:
            import pandas as pd
            import pyarrow as pa
        except ImportError:
            raise DeDuckError('Panda and PyArrow need to be imported in order to convert to dataframe')
        
        return pd.DataFrame.from_arrow(self._table)



def _get_gene_table(species):
    folder = resources.files('de_quack') / 'gene_tables'
    return folder / f'{species}_genes.parquet'

def _heal_genes(df, species):
    cursor = duckdb.connect(get_unique_conn())
    
    genes = pl.scan_parquet(_get_gene_table(species)).select('symbol', 'ensembl_id')
    df = df.lazy()
    df = df.with_columns(pl.col('gene_symbol').fill_null(pl.lit('')).alias('gene_symbol'))
    df = df.join(genes, on = 'ensembl_id', how = 'left')
    df = df.with_columns(pl.col('ensembl_id').fill_null(pl.lit('')).alias('ensembl_id'))
    df = df.join(
    genes,
    left_on="gene_symbol",  
    right_on="symbol",
    how="left", 
    suffix="_on_symbol"
    )
    final_df = df.select(pl.coalesce(['symbol', 'gene_symbol']).alias('gene_symbol'),
                                pl.coalesce(['ensembl_id', 'ensembl_id_on_symbol']).alias('ensembl_id'),
                                pl.all().exclude(['symbol', 'ensembl_id_on_symbol', 'gene_symbol', 'ensembl_id']))
    cursor.close()
    return final_df.collect()
        

class de_arrows:

    def __new__(cls, *args, metadata = None, ids=None):
        if metadata is None:
            return super().__new__(cls)
        if len(metadata) != len(args):
            if len(args) == 1 and isinstance(metadata, dict):
                return de_arrow(args[0], metadata=metadata)
        if len(args) == 1:
            return de_arrow(args[0], metadata=metadata[0])
        return super().__new__(cls)
    
    def __init__(self, *args, metadata = None, ids=None, heal_genes = False, species = None, keep_ids = False):
        """Initialize a de_arrows object with multiple differential expression tables.
        
        Parameters
        ----------
        *args : str or file paths
            Variable length argument list of table references (file paths or identifiers)
        metadata : list of dict
            List of metadata dictionaries, one for each table.
        """

        
        none_de_arrow = [table for table in args if not isinstance(table, de_arrow) and not isinstance(table, de_arrows)]
        if len(args) < 1:
            raise DeQuackError('At least one table must be provided')
        if metadata and not isinstance(metadata, list):
            if len(none_de_arrow) == 1:
                metadata = [metadata]
            else:
                raise DeQuackError('Metadata must be provided as a list of dictionaries')
        if len(metadata) != len(none_de_arrow):
            raise DeQuackError('Number of metadata dictionaries does not match number of tables')
        if keep_ids == True:
            if len(ids) != len(none_de_arrow):
                raise DeQuackError('Number of ids provided does not match number of tables')
            if isinstance(ids, int) and len(none_de_arrow) == 1:
                ids = [ids]
            if ids is None:
                ids = [i + 1 for i in range(len(args))]
        else:
            if ids is None:
                ids = [i + 1 for i in range(len(args))]
            if len(ids) != len(args):
                raise DeQuackError('Number of ids provided does not match number of tables')
        _check_ids(ids)
        
        
        self._table, self.experiment_metadata, ids = self.__class__._from_tables(*args, metadata=metadata, ids=ids, heal_genes = heal_genes, species = species)
        self.id = ids
        self.name = _get_experiment_attribute(self.experiment_metadata, 'experiment_name', ids)
        self.annotation_version = _get_experiment_attribute(self.experiment_metadata, 'annotation_version', ids)
        self.contrast = _get_experiment_attribute(self.experiment_metadata, 'contrast', ids)
        self.file = _get_experiment_attribute(self.experiment_metadata, 'file', ids)
        self.date = _get_experiment_attribute(self.experiment_metadata, 'date', ids)
        self.model = _get_experiment_attribute(self.experiment_metadata, 'model', ids)

    
    def __getattr__(self, name):
        # Delegate attribute access to the internal nanoarrow table without triggering recursion
        table = object.__getattribute__(self, '_table')
        return getattr(table, name)
    
    def __setattr__(self, name, value):
        if name in ('_table', 'experiment_metadata', 'name', 'annotation_version', 'contrast', 'file', 'date', 'id', 'model'):
            object.__setattr__(self, name, value)
        else:
            table = object.__getattribute__(self, '_table')
            setattr(table, name, value)
    
    def __getitem__(self, key):
        return self._table[key]
        



    
    def __len__(self):
        return len(self._table)
    
    def __repr__(self):
        return repr(self._table)
    
    def __str__(self):
        return(str(self._table))

    @classmethod
    def _from_tables(cls, *args, metadata = None, ids = None, heal_genes = False, species = None, keep_ids = False):
        experiment_columns=['experiment_id', 'model', 'date', 'file', 'experiment_name', 'contrast', 'annotation_version', 'normalization', 'extra_info']
        table = pl.DataFrame(schema = {'experiment_id': pl.Int32(), 'gene_symbol': pl.String(), 'ensembl_id': pl.String(), 'log2fc': pl.Float64(), 'logCPM': pl.Float64(), 'pvalue': pl.Float64(), 'padj': pl.Float64(), 'stat': pl.Float64(), 'other_info': pl.String()})
        flat = []
        final_ids = [id for id in ids]
        final_args = [arg for arg in args]
        if keep_ids == True:
            args = list(args)
            for arg in args:
                if isinstance(arg, de_arrow) or isinstance(arg, de_arrows):
                    table.extend(arg._table)
                    final_ids.insert(args.index(arg), arg.id)
                    for value in arg.experiment_metadata.values():
                        metadata.insert(args.index(arg), value)
                    args.remove(arg)
                    
        for arg, id in zip(args, ids):
            if isinstance(arg, de_arrow):
                arg = arg.replace_column(0, pl.lit(id).alias('experiment_id'))
                metadata.insert(args.index(arg), arg.experiment_metadata[arg.id])
            elif isinstance(arg, de_arrows):
                if not isinstance(id, list) and not isinstance(id, tuple) or len(id) != len(arg.id):
                    raise DeQuackError('ids must be a list or tuple for de_arrows objects and must be of same id length')
                id_lookup = pl.DataFrame({'old_id': arg.id, 'new_id': id}, schema = {'old_id': pl.Int32(), 'new_id': pl.Int32()})
                arg = arg.join(id_lookup, left_on = 'experiment_id', right_on = 'old_id', how = 'left')
                arg = arg.with_columns(pl.col('new_id').alias('experiment_id')).select(pl.all().exclude(['old_id', 'new_id']))
                for old_id in arg.id:
                    metadata.insert(args.index(arg), arg.experiment_metadata[old_id])
            else:
                arg = _to_polars_table(arg)
                arg = _order_columns(arg)
                arg = arg.insert_column(0, pl.lit(id).alias('experiment_id'))
            table.extend(arg)
        
        for id in final_ids:
            if isinstance(id, list) or isinstance(id, tuple):
                flat.extend(id)
            else:
                flat.append(id)
        _check_ids(flat)
        metadata_fields = _to_metadata(metadata, flat, final_args)
        return table, metadata_fields, flat
            
        

        
                
            
            
    @classmethod
    def _from_arrow(cls, table, metadata, ids, columns = None):
        instance = object.__new__(cls)
        instance._table = table
        instance.columns = _get_array_columns(table)
        instance.experiment_metadata = metadata
        instance.id = ids
        instance.name = _get_experiment_attribute(metadata, 'experiment_name', ids)
        instance.annotation_version = _get_experiment_attribute(metadata, 'annotation_version', ids)
        instance.contrast = _get_experiment_attribute(metadata, 'contrast', ids)
        instance.file = _get_experiment_attribute(metadata, 'file', ids)
        instance.date = _get_experiment_attribute(metadata, 'date', ids)
        instance.model = _get_experiment_attribute(metadata, 'model', ids)
        if columns is not None:
            instance.columns = columns
        else:
            instance.columns = _get_array_columns(table)
        return instance

    def _set_id_dict(self, ids: dict):
        for old_id, new_id in ids.items():
            if old_id not in self.id:
                raise DeQuackError(f'{old_id} is not in existing ids. Existing ids are {self.id}')
        _check_ids(list(ids.values()))
        de = de_quackling(get_unique_conn()).connect()
        try:
            old_ids = list(ids.keys())
            new_ids = list(ids.values())
            ids_rel = de.conn.sql(f'SELECT UNNEST({old_ids}) AS old_id, UNNEST({new_ids}) AS new_id')
            arrow_table = self._table
            rel = de.conn.sql(_de_arrows_queries['set_ids'])
            arrow_table = na.Array(rel)
            updated_metadata = {new_id : self.experiment_metadata.get(old_id, {}) for old_id, new_id in ids.items()}
            for old_id in self.experiment_metadata.keys():
                if old_id not in updated_metadata.keys() and old_id not in old_ids:
                    updated_metadata[old_id] = self.experiment_metadata[old_id]
            arrow_table = self._from_arrow(arrow_table, updated_metadata, ids.values())
            return arrow_table
        finally:
            de.conn.close()
        

    def _set_id_list(self, ids: list):
        print(self.id)
        if len(ids) != len(self.id):
            raise DeQuackError('Number of provided ids does not match number of existing ids')
        _check_ids(ids)
            
        de=de_quackling(get_unique_conn()).connect()
        try:
            ids_rel = de.conn.sql(f'SELECT UNNEST({self.id}) AS old_id, UNNEST({ids}) AS new_id')
            arrow_table = self._table
            rel = de.conn.sql(_de_arrows_queries['set_ids'])
            arrow_table = na.Array(rel)
            updated_metadata = {new_id: self.experiment_metadata.get(old_id, {}) for old_id, new_id in zip(self.id, ids)}
            arrow_table = self._from_arrow(arrow_table, updated_metadata, ids)
            return arrow_table
        finally:
            de.conn.close()
    
    def set_id(self, ids):
        if isinstance(ids, dict):
            return self._set_id_dict(ids)
        if isinstance(ids, list):
            return self._set_id_list(ids)
        raise DeQuackError(f'set_ids expected a dictionary or a list, not {type(ids)}')
    
    def get_upregulated(self, log2fc=0, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        de=de_quackling(get_unique_conn()).fast_connect()
        try:
            arrow_table=self._table
            rel=de.conn.sql(_de_arrow_queries['get_upregulated'], params={'log2fc':log2fc, 'pvalue':pvalue, 'logCPM': logCPM})
            arrow_table=na.Array(rel)
            arrow_table = self._finalize_table(arrow_table)
            return arrow_table
        finally:
            de.conn.close()

    def get_significant_genes(self, log2fc=1, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        de=de_quackling(get_unique_conn()).fast_connect()
        try:
            arrow_table=self._table
            rel=de.conn.sql(_de_arrow_queries['get_important_genes'], params={'log2fc':log2fc, 'pvalue':pvalue, 'logCPM': logCPM})
            arrow_table=na.Array(rel)
            arrow_table=self._finalize_table(arrow_table)
            return arrow_table
        finally:
            de.conn.close()
    
    def get_downregulated(self, log2fc=0, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        de=de_quackling(get_unique_conn()).fast_connect()
        try:
            arrow_table=self._table
            rel=de.conn.sql(_de_arrow_queries['get_downregulated'], params={'log2fc':log2fc, 'pvalue':pvalue, 'logCPM': logCPM})
            arrow_table = na.Array(rel)
            arrow_table = self._finalize_table(arrow_table)
            return arrow_table
        finally:
            de.conn.close()

    def _finalize_table(self, arrow_table):
        de = de_quackling(get_unique_conn()).fast_connect()
        try:
            ids = de.conn.sql("SELECT DISTINCT experiment_id FROM arrow_table").fetchall()
            ids = [id[0] for id in ids]
            new_metadata = {}
            for id in ids:
                new_metadata[id] = self.experiment_metadata.get(id, '')
            if len(ids) < 2:
                if ids is None:
                    return de_arrow._from_arrow(arrow_table, {}, None)
                return de_arrow._from_arrow(arrow_table, self.experiment_metadata[ids[0]], id)
            return self._from_arrow(arrow_table, new_metadata, ids)
        finally:
            de.conn.close()
    
    def add_experiment(self, data, metadata = None, id = None):
        if isinstance(data, de_arrow):
            return self._add_experiment_arrow(data, metadata, id)
        if isinstance(data, de_arrows):
            if metadata is not None:
                if not isinstance(metadata, list) or not isinstance(metadata[0], dict):
                    raise DeQuackError('Metadata must be a list of dictionaries for de_arrows objects!')
                if len(id) != len(data.id):
                    raise DeQuackError('Length of id list does not match length of incoming de_arrows object')
            return self._add_experiment_arrows(data, metadata, id)
        return self._add_experiment_data(data, metadata, id)

    
    def _add_experiment_arrows(self, data, metadata = None, ids = None):
        de = de_quackling(get_unique_conn()).fast_connect()
        try:
            arrow_table = self._table
            if ids is None:
                ids = data.id
                _check_ids(ids + self.id)
            else:
                _check_ids(ids + self.id)
            self_metadata = self.experiment_metadata
            if metadata is None:
                new_metadata = data.experiment_metadata
                for id, value in zip(ids, new_metadata.values()):
                    self_metadata[id] = value
            else:
                new_metadata = metadata
                for meta, id in zip(ids, new_metadata):
                    m, extra_info = ExperimentMetadata.to_dict(new_metadata)
                    for key, value in extra_info.items():
                        m[key] = value
                    self_metadata[n] = m
            arrow_table = de.conn.sql("SELECT * FROM arrow_table UNION ALL SELECT * FROM data")
            arrow_table = na.Array(arrow_table)
            ids = [id for id in self_metadata.keys()]
            return self._from_arrow(arrow_table, self_metadata, ids)
        finally:
            de.conn.close()
        
        
    def _add_experiment_arrow(self, data, metadata = None, id = None):
        arrow_table = self._table
        if not isinstance(id, int) and id is not None:
            raise DeQuackError('id has to be an integer')
        de = de_quackling(get_unique_conn()).fast_connect()
        try:
            new_metadata = ExperimentMetadata.to_dict(metadata) if metadata is not None else data.experiment_metadata
            id = id if id is not None else data.id
            if id in self.id:
                raise DeQuackError('Duplicate ids found')
            metadata = self.experiment_metadata
            arrow_table = de.conn.sql("SELECT * FROM arrow_table UNION ALL SELECT * FROM data")
            arrow_table = na.Array(arrow_table)
            metadata[id] = list(new_metadata.values())[0]
            ids = self.id
            ids.append(id)
            return self._from_arrow(arrow_table, metadata, ids)
        finally:
            de.conn.close()

    def _add_experiment_data(self, data, metadata = None, id = None):
        arrow_table = self._table
        if metadata is None:
            raise DeQuackError('No metadata provided for the data')
        if id is None:
            raise DeQuackError('No id provided')
        de = de_quackling(get_unique_conn()).connect()
        try:
            new_metadata, extra_info=ExperimentMetadata().to_dict(metadata, data)
            for key, value in json.loads(extra_info).items():
                new_metadata[key]=value
            new_metadata={id: new_metadata}
            de._preprocess(data)
            de_arrow_insertion_view = de._create_temp_view()
            data = de.conn.sql(_de_arrow_queries['insert_to_de_arrow'], params={'id':id})
            if id in self.id:
                raise DeQuackError('Duplicate ids found')
            metadata = self.experiment_metadata
            arrow_table = de.conn.sql("SELECT * FROM arrow_table UNION ALL SELECT * FROM data")
            arrow_table = na.Array(arrow_table)
            metadata[id] = new_metadata[id]
            ids = self.id
            ids.append(id)
            return self._from_arrow(arrow_table, metadata, ids)
        finally:
            de.conn.close()
    
    def get_experiment(self, id = None, name = None, model = None, annotation_version = None,  normalization = None, date = None, contrast = None, file = None):
        try:
            date = datetime.datetime.strptime(date, '%Y-%m-%d').date() if date else None
        except ValueError:
            raise ProcessingError(f"Invalid date format: {date}. Expected format is YYYY-MM-DD.")
        metadata = self.experiment_metadata
        selected_metadata = set(self.id)
        not_selected = set([])
        if id is not None:
            if isinstance(id, list):
                selected_metadata = set([i for i in self.id if i in id])
            elif isinstance(id, int):
                selected_metadata = {id}
        if name is not None:
            if isinstance(name, list):
                name = [n.lower().strip() for n in name]
                not_selected.update([i for i, n in self.name.items() if n.lower.strip() not in name])
            elif isinstance(name, str):
                not_selected.update([i for i, n in self.name.items() if n.lower() not in name.lower()])
        if model is not None:
            if isinstance(model, list):
                model = [n.lower().strip() for n in model]
                not_selected.update([i for i, n in self.model.items() if n.lower.strip() not in model])
            elif isinstance(model, str):
                not_selected.update([i for i, n in self.model.items() if n.lower() not in model.lower()])
        if annotation_version is not None:
            if isinstance(annotation_version, list):
                annotation_version = [n.lower().strip() for n in annotation_version]
                not_selected.update([i for i, n in self.annotation_version.items() if n.lower.strip() not in annotation_version])
            elif isinstance(annotation_version, str):
                not_selected.update([i for i, n in self.annotation_version.items() if n.lower() not in annotation_version.lower()])
        if normalization is not None:
            if isinstance(normalization, list):
                normalization = [n.lower().strip() for n in normalization]
                not_selected.update([i for i, n in self.normalization.items() if n.lower.strip() not in normalization])
            elif isinstance(normalization, str):
                not_selected.update([i for i, n in self.normalization.items() if n.lower() not in normalization.lower()])
        if date is not None:
            if isinstance(date, list):
                selected_metadata.update([i for i, n in self.date.items() if n not in date])
            elif isinstance(date, datetime.date):
                not_selected.update([i for i, n in self.date.items() if n != date])
        if contrast is not None:
            if isinstance(contrast, list):
                contrast = [n.lower().strip() for n in contrast]
                not_selected.update([i for i, n in self.contrast.items() if n.lower.strip() not in contrast])
            elif isinstance(contrast, str):
                not_selected.update([i for i, n in self.contrast.items() if n.lower() not in contrast.lower()])
        if file is not None:
            if isinstance(file, list):
                file = [n.lower().strip() for n in file]
                not_selected.update([i for i, n in self.file.items() if n.lower.strip() not in file])
            elif isinstance(file, str):
                not_selected.update([i for i, n in self.file.items() if n.lower() not in file.lower()])
        
        selected_metadata = list(selected_metadata - not_selected)
        if len(selected_metadata) == 0:
            raise DeQuackError('No experiments found with the provided metadata')

        de = de_quackling(get_unique_conn()).fast_connect()
        try:
            arrow_table = self._table
            rel = de.conn.sql(_de_arrows_queries['get_experiment'], params = {
                'ids': selected_metadata
            })
            arrow_table = na.Array(rel)
            return self._finalize_table(arrow_table)
        finally:
            de.conn.close()

    def get_gene(self, gene_symbol = None, ensembl_id = None, id = None):
        required_columns = {'gene_symbol', 'ensembl_id'}
        _check_columns(required_columns, self.columns)
        de = de_quackling(get_unique_conn()).fast_connect()
        try:
            arrow_table = self._table
            rel = de.conn.sql(_de_arrows_queries['get_gene'], 
            params = {
            'gene_symbol': gene_symbol,
            'ensembl_id': ensembl_id,
            'id': id
            })
            arrow_table = na.Array(rel)
            return self._finalize_table(arrow_table)
        finally:
            de.conn.close()
    
    def insert(self, file):
        required_columns = {'padj', 'pvalue', 'log2fc', 'gene_symbol', 'ensembl_id', 'logCPM', 'stat', 'other_info', 'experiment_id'}
        _check_columns(required_columns, self.columns)
        if not isinstance(file, str):
            raise TypeError(f'file must be a string, not {type(file)}')
        arrow_table = self._table
        metadata = self.experiment_metadata
        new_meta = []
        for meta in metadata.values():
            m, other_info = ExperimentMetadata().to_dict(meta, None)
            new_meta.append((m, other_info))
        de = de_quackling(file).connect()
        try:
            tables = []
            for id in self.id:
                rel = de.conn.sql('SELECT * FROM arrow_table WHERE experiment_id = $id', params = {'id': id})
                tables.append(na.Array(rel))
            for meta, table in zip(new_meta, tables):
                arrow = table
                sample_data = de.conn.execute('SELECT * FROM (SELECT * FROM arrow LIMIT 100)').fetchall()
                data_sig = hashlib.md5(str(sample_data).encode()).hexdigest()
                id = de._create_experiment(meta[0], data_sig, meta[1])
                de.conn.sql('INSERT INTO gene_results (experiment_id, gene_symbol, ensembl_id, log2fc, logCPM, pvalue, padj, stat, other_info) SELECT * REPLACE($id AS experiment_id) FROM arrow', params = {'id': id})
            de.conn.commit()
            return self
        finally:
            de.conn.close()

            


        
    
    def df(self):
        try:
            import pandas as pd
            import pyarrow as pa
        except ImportError:
            raise DeDuckError('Panda and PyArrow need to be imported in order to convert to dataframe')
        
        return pd.DataFrame.from_arrow(self._table)

    


def get_unique_conn():
    unique_id = uuid.uuid4().hex[:8]
    return(f':memory:de_quack_{unique_id}')

def _get_experiment_attribute(metadata, field, ids):
    if isinstance(ids, list):
        return {id: metadata.get(id, {}).get(field) for id in ids}
    else:
        return {ids: metadata.get(ids, {}).get(field)}

def _check_columns(required_columns, columns):
        missing_columns = required_columns - set(columns)
        if missing_columns:
            raise DeQuackError(f'De_arrow object is missing {missing_columns} for this function')

def _check_ids(ids):
    for id in ids:
        if ids.count(id) > 1:
            raise DeQuackError(f'id {id} occurs multiple times in object ids')
        if isinstance(id, str):
            raise TypeError('String found in id list')

def _to_metadata(metadata, id, tables):
    meta_list = []
    to_remove_tables = []
    to_remove_metadata = []
    print(metadata)
    for table in tables:
        if isinstance(table, de_arrows):
            to_remove_tables.append(table)
            print(table.experiment_metadata)
            for n in table.id:
                print(metadata.index(table.experiment_metadata[n]) + table.id.index(n))
                m = metadata[metadata.index(table.experiment_metadata[n]) + table.id.index(n)]
                meta_list.append(m)
                to_remove_metadata.append(m)
    for table in to_remove_tables:
        tables.remove(table)
    for meta in to_remove_metadata:
        metadata.remove(meta)
    for meta, table in zip(metadata, tables):
        print(meta, table)
        if isinstance(table, de_arrow):
            meta_list.append(meta)
        else:
            m, extra_info = ExperimentMetadata().to_dict(meta, table)
            extra_info = json.loads(extra_info)
            for key, value in extra_info.items():
                m[key] = value
            meta_list.append(m)
    
    metadata_dict = {}
    for meta, id in zip(meta_list, id):
        metadata_dict[id] = meta
    return metadata_dict

            

def _to_polars_table(table):
    if isinstance(table, str):
        path = os.path.abspath(table)
        if os.path.exists(path):
            if path.endswith('.parquet'):
                return pl.read_parquet(path)
            else:
                with open(path, 'r') as f:
                    first_line = f.readline()
                    separator = "\t" if "\t" in first_line else (";" if ";" in first_line else ",")
                return pl.read_csv(path, separator=separator)
    class_name = table.__class__.__name__
    if class_name == 'Table':
        return pl.from_arrow(table)
    if class_name == 'DataFrame':
        return pl.from_pandas(table)
    raise ValueError(f"Unsupported table type: {class_name}. Supported types are: str (file path), pyarrow.Table, pandas.DataFrame.")

def _order_columns(df):
    columns = df.columns
    extra_columns = [col for col in columns if col not in _GENE_ALIAS_TO_COLUMN.values()]
    real_columns = [_GENE_ALIAS_TO_COLUMN.get(col.lower()) for col in columns if _GENE_ALIAS_TO_COLUMN.get(col.lower())]
    essential_cols = {'padj', 'log2fc', 'padj'}
    if essential_cols - set(real_columns):
        raise DeQuackError(f'Missing essential columns: {essential_cols - set(real_columns)}')
    df = df.rename({col: _GENE_ALIAS_TO_COLUMN.get(col.lower()) for col in columns if _GENE_ALIAS_TO_COLUMN.get(col.lower())})
    if 'gene_name' in df.columns:
        sample = df.select('gene_name').to_series().head(1).to_list()
        if re.match('^ENS', sample[0]):
            if 'ensembl_id' not in df.columns:
                df = df.rename({'gene_name': 'ensembl_id'})                 
        else:
            if 'gene_symbol' not in df.columns:
                df = df.rename({'gene_name': 'gene_symbol'})
    expressions = []
    for col in _gene_columns:
        if col in df.columns:
            expressions.append(pl.col(col))
        else:
            expressions.append(pl.lit(None).alias(col))
    extra_columns = set([col for col in df.columns if col not in _gene_columns and col != 'experiment_id'])
    if not extra_columns:
        expressions.append(pl.lit(None).alias('other_info'))
    else:
        expressions.append(pl.struct(extra_columns).struct.json_encode().alias('other_info'))
    df = df.select(expressions)

    return df



