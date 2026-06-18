import nanoarrow as na
from .core import de_quackling
from .utilities import ExperimentMetadata, DE_ARROW_QUERIES, DE_ARROWS_QUERIES
import duckdb
from .exceptions import DeQuackError, ProcessingError, DuplicateExperimentError, DuplicateGeneTableError
import json
import uuid
import re
_de_arrow_queries=DE_ARROW_QUERIES
_de_arrows_queries=DE_ARROWS_QUERIES
_gene_columns = ['padj', 'p_value', 'log2fc', 'gene_symbol', 'ensembl_id', 'logCPM', 'experiment_id']
class de_arrow:
    
    def __init__(self, info, id = None, metadata=None, **fields):
        
        if metadata is None:
            metadata = fields
        if len(metadata) == 0:
            raise ProcessingError('No metadata provided in de_arrow object')
        
        self._table, self.experiment_metadata = self.__class__._to_de_arrow(info, metadata, id)
        self.columns = _get_array_columns(self._table)
        self.name = _get_experiment_attribute(self.experiment_metadata, 'name', id)
        self.annotation_version = _get_experiment_attribute(self.experiment_metadata, 'annotation_version', id)
        self.contrast = _get_experiment_attribute(self.experiment_metadata, 'contrast', id)
        self.file = _get_experiment_attribute(self.experiment_metadata, 'file', id)
        self.date = _get_experiment_attribute(self.experiment_metadata, 'date', id)
        self.id = id

        
    
    def __getattr__(self, name):
        # Delegate attribute access to the internal nanoarrow table without triggering recursion
        table = object.__getattribute__(self, '_table')
        return getattr(table, name)
    
    def __setattr__(self, name, value):
        if name in ('_table', 'experiment_metadata', 'name', 'annotation_version', 'contrast', 'file', 'date', 'id'):
            object.__setattr__(self, name, value)
        else:
            table = object.__getattribute__(self, '_table')
            setattr(table, name, value)
    
    def __getitem__(self, key):
        if isinstance(key, int):
            return self._table[key]
        elif isinstance(key, str):
            if key in _gene_columns:
                de = de_quackling(get_unique_conn()).fast_connect()
                arrow_table=self._table
                arrow_table = de.conn.sql(f"SELECT {key} FROM arrow_table")
                arrow_table = na.Array(arrow_table)
                arrow_table = self._from_arrow(arrow_table, self.experiment_metadata, self.id)
                return arrow_table
            else:
                raise KeyError(f'{key} is not in the table')
        elif isinstance(key, list):
            for column in key:
                if column not in _gene_columns:
                    raise KeyError(f'{column} is not in the table')
            de = de_quackling(get_unique_conn()).fast_connect()
            arrow_table = self._table
            arrow_table = de.conn.sql(f"SELECT COLUMNS({key}) FROM arrow_table")
            arrow_table = na.Array(arrow_table)
            arrow_table = self._from_arrow(arrow_table, self.experiment_metadata, self.id)
            return arrow_table

    def __len__(self):
        return len(self._table)
    
    def __repr__(self):
        return repr(self._table)
    
    def __str__(self):
        return str(self._table)
    

    
    @classmethod
    def _to_de_arrow(self, info, metadata, experiment_id=None, heal_genes=False):
        de=de_quackling(get_unique_conn()).connect()
        metadata_fields, extra_info=ExperimentMetadata().to_dict(metadata, info)
        for key, value in json.loads(extra_info).items():
            metadata_fields[key]=value
        metadata_fields={experiment_id: metadata_fields}
        de._preprocess(info)
        de_arrow_insertion_view=de._create_temp_view()
        
        rel=de.conn.sql(_de_arrow_queries['insert_to_de_arrow'], params={'id':experiment_id})
        schema=na.Schema(na.Type.INT64, metadata={f'{key}.exp{experiment_id}': value for key, value in metadata_fields[experiment_id].items()})
        table=na.Array(rel, schema)
        
        
        de.conn.close()
        return table, metadata_fields
    
    @classmethod
    def _from_arrow(cls, arrow_array, metadata, id):
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
        instance.columns = _get_array_columns(arrow_array)
        instance.experiment_metadata = metadata
        instance.name = _get_experiment_attribute(metadata, 'name', id)
        instance.annotation_version = _get_experiment_attribute(metadata, 'annotation_version', id)
        instance.contrast = _get_experiment_attribute(metadata, 'contrast', id)
        instance.file = _get_experiment_attribute(metadata, 'file', id)
        instance.date = _get_experiment_attribute(metadata, 'date', id)
        instance.id = id
        return instance
    
    
    def get_significant_genes(self, log2fc=1, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        de=de_quackling(get_unique_conn()).fast_connect()
        arrow_table=self._table
        rel=de.conn.sql(_de_arrow_queries['get_important_genes'], params={'log2fc':log2fc, 'pvalue':pvalue, 'logCPM': logCPM})
        arrow_table=na.Array(rel)
        arrow_table=self._from_arrow(arrow_table, self.experiment_metadata, self.id)
        de.conn.close()
        return arrow_table
    
    def get_downregulated(self, log2fc=0, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        de=de_quackling(get_unique_conn()).fast_connect()
        arrow_table=self._table
        rel=de.conn.sql(_de_arrow_queries['get_downregulated'], params={'log2fc':log2fc, 'pvalue':pvalue, 'logCPM': logCPM})
        arrow_table=self._from_arrow(arrow_table, self.experiment_metadata, self.id)
        de.conn.close()
        return arrow_table
    
    def get_upregulated(self, log2fc=0, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        de=de_quackling(get_unique_conn()).fast_connect()
        arrow_table=self._table
        rel=de.conn.sql(_de_arrow_queries['get_upregulated'], params={'log2fc':log2fc, 'pvalue':pvalue, 'logCPM': logCPM})
        arrow_table=na.Array(rel)
        arrow_table=self._from_arrow(arrow_table, self.experiment_metadata, self.id)
        de.conn.close()
        return arrow_table
    
    def set_id(self, id):
        de=de_quackling(get_unique_conn()).fast_connect()
        arrow_table=self._table
        rel=de.conn.sql(_de_arrow_queries['set_arrow_id'], params={'id':id})
        arrow_table=na.Array(rel)
        old_id = self.id
        updated_metadata = {id: self.experiment_metadata.get(old_id, {})}
        arrow_table=self._from_arrow(arrow_table, updated_metadata, id)
        de.conn.close()
        return arrow_table

    
    def df(self):
        try:
            import pandas as pd
            import pyarrow as pa
        except ImportError:
            raise DeDuckError('Panda and PyArrow need to be imported in order to convert to dataframe')
        
        return pd.DataFrame.from_arrow(self._table)
        

class de_arrows:

    def __new__(cls, *args, metadata, ids=None):
        if len(metadata) != len(args):
            if len(args) == 1 and isinstance(metadata, dict):
                return de_arrow(args[0], metadata=metadata)
            else:
                raise DeQuackError(f'Metadata length of {len(metadata)} does not match number of tables ({len(args)})')
        if len(args) == 1:
            return de_arrow(args[0], metadata=metadata[0])
        return super().__new__(cls)
    
    def __init__(self, *args, metadata, ids=None):
        """Initialize a de_arrows object with multiple differential expression tables.
        
        Parameters
        ----------
        *args : str or file paths
            Variable length argument list of table references (file paths or identifiers)
        metadata : list of dict
            List of metadata dictionaries, one for each table.
        """
        if not isinstance(metadata, list):
            raise DeQuackError('Metadata must be provided as a list of dictionaries')
        if len(args) < 1:
            raise DeQuackError('At least one table must be provided')
        if not metadata:
            raise DeQuackError('Metadata must be provided as a list of dictionaries')
        if ids is not None:
            if len(ids) != len(args):
                raise DeQuackError('Number of ids provided does not match number of tables')
            for id in ids:
                if isinstance(id, str):
                    raise DeQuackError('type str was found in id list')
        else:
            ids = [i for i in range(len(args))]
        
        self._table, self.experiment_metadata = self.__class__._from_tables(*args, metadata=metadata, ids=ids)
        self.id = ids
        self.name = _get_experiment_attribute(self.experiment_metadata, 'name', ids)
        self.columns = _get_array_columns(self._table)
        self.annotation_version = _get_experiment_attribute(self.experiment_metadata, 'annotation_version', ids)
        self.contrast = _get_experiment_attribute(self.experiment_metadata, 'contrast', ids)
        self.file = _get_experiment_attribute(self.experiment_metadata, 'file', ids)
        self.date = _get_experiment_attribute(self.experiment_metadata, 'date', ids)

    
    def __getattr__(self, name):
        # Delegate attribute access to the internal nanoarrow table without triggering recursion
        table = object.__getattribute__(self, '_table')
        return getattr(table, name)
    
    def __setattr__(self, name, value):
        if name in ('_table', 'experiment_metadata', 'name', 'annotation_version', 'contrast', 'file', 'date', 'id'):
            object.__setattr__(self, name, value)
        else:
            table = object.__getattribute__(self, '_table')
            setattr(table, name, value)
    
    def __getitem__(self, key):
        if isinstance(key, int):
            return self._table[key]
        elif isinstance(key, str):
            if key in _gene_columns:
                de = de_quackling(get_unique_conn()).fast_connect()
                arrow_table=self._table
                arrow_table = de.conn.sql(f"SELECT {key} FROM arrow_table")
                arrow_table = na.Array(arrow_table)
                arrow_table = self._from_arrow(arrow_table, self.experiment_metadata, self.id)
                return arrow_table
            else:
                raise KeyError(f'{key} is not in the table')
        elif isinstance(key, list):
            for column in key:
                if column not in _gene_columns:
                    raise KeyError(f'{column} is not in the table')
            de = de_quackling(get_unique_conn()).fast_connect()
            arrow_table = self._table
            arrow_table = de.conn.sql(f"SELECT COLUMNS({key}) FROM arrow_table")
            arrow_table = na.Array(arrow_table)
            arrow_table = self._from_arrow(arrow_table, self.experiment_metadata, self.id)
            return arrow_table



    
    def __len__(self):
        return len(self._table)
    
    def __repr__(self):
        return repr(self._table)
    
    def __str__(self):
        return str(self._table)

    @classmethod
    def _from_tables(cls, *args, metadata, ids = None):
        experiment_columns=['experiment_id', 'model', 'date', 'file', 'experiment_name', 'contrast', 'annotation_version', 'normalization', 'extra_info']
        
        de = de_quackling(get_unique_conn()).connect()
        try:
            for table, meta in zip(args, metadata):
                de.ingest(table, meta)
    
            ids=de.conn.sql(f'SELECT UNNEST({ids}) AS ids')
            rel = de.conn.sql(_de_arrows_queries['insertion_from_table'])
            
            
            
            meta_rel = rel.select(', '.join(experiment_columns)).distinct()
            rel = rel.select('experiment_id, gene_symbol, ensembl_id, log2fc, logCPM, pvalue, padj, other_info')
            metadata_list = [dict(zip(experiment_columns, row)) for row in meta_rel.fetchall()]
            
            # Parse extra_info JSON and flatten into metadata
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
            
            schema = na.Schema(na.Type.INT64)
            table = na.Array(rel, schema)
            
           
            return table, metadata_fields
            
        finally:
            de.conn.close()
    
    @classmethod
    def _from_arrow(cls, table, metadata, ids):
        instance = object.__new__(cls)
        instance._table = table
        instance.columns = _get_array_columns(table)
        instance.experiment_metadata = metadata
        instance.id = ids
        instance.name = _get_experiment_attribute(metadata, 'name', ids)
        instance.annotation_version = _get_experiment_attribute(metadata, 'annotation_version', ids)
        instance.contrast = _get_experiment_attribute(metadata, 'contrast', ids)
        instance.file = _get_experiment_attribute(metadata, 'file', ids)
        instance.date = _get_experiment_attribute(metadata, 'date', ids)
        return instance

    def set_ids(self, ids: list):
        if not isinstance(ids, list):
            raise DeQuackError('ids must be provided in a list')
        if len(ids) != len(self.id):
            raise DeQuackError('Number of provided ids does not match number of existing ids')
        for id in ids:
            if isinstance(id, str):
                raise DeQuackError('String found in id')
            
        de=de_quackling(get_unique_conn()).connect()
        ids_rel = de.conn.sql(f'SELECT UNNEST({self.id}) AS old_id, UNNEST({ids}) AS new_id')
        arrow_table = self._table
        rel = de.conn.sql(_de_arrows_queries['set_ids'])
        arrow_table = na.Array(rel)
        updated_metadata = {new_id: self.experiment_metadata.get(old_id, {}) for old_id, new_id in zip(self.id, ids)}
        arrow_table = self._from_arrow(arrow_table, updated_metadata, ids)
        return arrow_table
    
    def get_upregulated(self, log2fc=0, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        de=de_quackling(get_unique_conn()).fast_connect()
        arrow_table=self._table
        rel=de.conn.sql(_de_arrow_queries['get_upregulated'], params={'log2fc':log2fc, 'pvalue':pvalue, 'logCPM': logCPM})
        arrow_table=na.Array(rel)
        arrow_table = self._finalize_table(arrow_table)
        de.conn.close()
        return arrow_table

    def get_significant_genes(self, log2fc=1, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        de=de_quackling(get_unique_conn()).fast_connect()
        arrow_table=self._table
        rel=de.conn.sql(_de_arrow_queries['get_important_genes'], params={'log2fc':log2fc, 'pvalue':pvalue, 'logCPM': logCPM})
        arrow_table=na.Array(rel)
        arrow_table=self._finalize_table(arrow_table)
        de.conn.close()
        return arrow_table
    
    def get_downregulated(self, log2fc=0, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        de=de_quackling(get_unique_conn()).fast_connect()
        arrow_table=self._table
        rel=de.conn.sql(_de_arrow_queries['get_downregulated'], params={'log2fc':log2fc, 'pvalue':pvalue, 'logCPM': logCPM})
        arrow_table = na.Array(rel)
        arrow_table = self._finalize_table(arrow_table)
        de.conn.close()
        return arrow_table

    def _finalize_table(self, arrow_table):
        de = de_quackling(get_unique_conn()).fast_connect()
        ids = de.conn.sql("SELECT COUNT(DISTINCT experiment_id) FROM arrow_table").fetchone()
        if ids[0] < 2:
            id = de.conn.sql("SELECT DISTINCT experiment_id FROM arrow_table").fetchone()
            de.conn.close()
            return de_arrow._from_arrow(arrow_table, self.experiment_metadata[id[0]], id)
        de.conn.close()
        return self._from_arrow(arrow_table, self.experiment_metadata, self.id)
        
        


        
    
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

def _get_array_columns(table):
    schema = table.schema
    return [field.name for field in schema.fields]

def _check_columns(required_columns, columns):
        missing_columns = required_columns - set(columns)
        if missing_columns:
            raise DeQuackError(f'De_arrow object is missing {missing_columns} for this function')

            







