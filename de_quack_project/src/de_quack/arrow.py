import json
import os
import re
import time
import uuid
from importlib import resources
import hashlib
import duckdb
import polars as pl
import polars.selectors as cs
from .core import de_quackling
from .exceptions import DeQuackError, ProcessingError, DuplicateGeneTableError
from .utilities import DE_ARROW_QUERIES, DE_ARROWS_QUERIES, ExperimentMetadata, gene_columns, gene_mapping, CORE_QUERIES


_de_arrow_queries = DE_ARROW_QUERIES
_de_arrows_queries = DE_ARROWS_QUERIES
_core_queries = CORE_QUERIES
_gene_mapping_queries = gene_mapping
_GENE_ALIAS_TO_COLUMN = {
    alias.lower(): canonical
    for canonical, aliases in gene_columns.items()
    for alias in aliases + [canonical]
}
_gene_columns = ['gene_symbol', 'ensembl_id', 'log2fc', 'logCPM', 'pvalue', 'padj', 'stat']


def get_unique_conn():
    unique_id = uuid.uuid4().hex[:8]
    return f':memory:de_quack_{unique_id}'

def _clean_df(df):
    df = df.with_columns(pl.col('gene_symbol').fill_null(pl.lit('')).alias('gene_symbol'),
                        pl.col('ensembl_id').fill_null(pl.lit('')).alias('ensembl_id'),
                        pl.col('log2fc').fill_null(pl.lit(0.0)).alias('log2fc'),
                        pl.col('pvalue').fill_null(pl.lit(1.0)).alias('pvalue'),
                        pl.col('logCPM').fill_null(pl.lit(0.0)).alias('logCPM'),
                        pl.col('other_info').fill_null(pl.lit('{}')).alias('other_info'),
                        pl.col('padj').fill_null(pl.lit(1.0)).alias('padj'))
    return df


def _to_polars_table(table):
    if isinstance(table, pl.DataFrame):
        return table.clone()
    if isinstance(table, pl.LazyFrame):
        return table.collect()
    if hasattr(table, 'to_frame') and callable(table.to_frame):
        return _to_polars_table(table.to_frame())

    try:
        import pandas as pd
    except ImportError:
        pd = None

    if pd is not None and isinstance(table, pd.DataFrame):
        return pl.from_pandas(table)

    if isinstance(table, str):
        path = os.path.abspath(table)
        if not os.path.exists(path):
            raise FileNotFoundError(f'File not found: {path}')
        if path.lower().endswith('.parquet'):
            return pl.read_parquet(path)
        with open(path, 'r', encoding='utf-8') as handle:
            first_line = handle.readline()
        separator = '\t' if '\t' in first_line else (';' if ';' in first_line else ',')
        return pl.read_csv(path, separator=separator)

    if isinstance(table, dict):
        return pl.DataFrame(table)

    if isinstance(table, list):
        return pl.DataFrame(table)

    return pl.DataFrame(table)


def _check_columns(required_columns, columns):
    missing_columns = required_columns - set(columns)
    if missing_columns:
        raise DeQuackError(f'De_arrow object is missing {missing_columns} for this function')


def _check_ids(ids):
    for item in ids:
        if ids.count(item) > 1:
            raise DeQuackError(f'id {item} occurs multiple times in object ids')
        if isinstance(item, str):
            raise TypeError('String found in id list')


def _get_experiment_attribute(metadata, field, ids):
    if isinstance(ids, list):
        return {item: metadata.get(item, {}).get(field) for item in ids}
    return {ids: metadata.get(ids, {}).get(field)}


def _get_gene_table(species):
    folder = resources.files('de_quack') / 'gene_tables'
    return folder / f'{species}_genes.parquet'


def _heal_genes(df, species):
    df = _to_polars_table(df).lazy()
    genes = pl.scan_parquet(_get_gene_table(species)).select('symbol', 'ensembl_id')

    df = df.with_columns(pl.col('gene_symbol').fill_null(pl.lit('')).alias('gene_symbol'))
    df = df.join(genes, on='ensembl_id', how='left')
    df = df.with_columns(pl.col('ensembl_id').fill_null(pl.lit('')).alias('ensembl_id'))
    df = df.join(genes, left_on='gene_symbol', right_on='symbol', how='left', suffix='_on_symbol')

    return df.select(
        pl.coalesce(['symbol', 'gene_symbol']).alias('gene_symbol'),
        pl.coalesce(['ensembl_id', 'ensembl_id_on_symbol']).alias('ensembl_id'),
        pl.all().exclude(['symbol', 'ensembl_id_on_symbol', 'gene_symbol', 'ensembl_id'])
    ).collect()


def _order_columns(df):
    df = _to_polars_table(df)
    rename_map = {
        column: _GENE_ALIAS_TO_COLUMN[column.lower()]
        for column in df.columns
        if column.lower() in _GENE_ALIAS_TO_COLUMN and column != 'experiment_id'
    }
    if rename_map:
        df = df.rename(rename_map)

    if 'gene_name' in df.columns:
        sample = df.select('gene_name').to_series().head(1).to_list()
        if sample and sample[0] is not None and re.match(r'^ENS', str(sample[0])):
            if 'ensembl_id' not in df.columns:
                df = df.rename({'gene_name': 'ensembl_id'})
        elif 'gene_symbol' not in df.columns:
            df = df.rename({'gene_name': 'gene_symbol'})

    expressions = []
    if 'experiment_id' in df.columns:
        expressions.append(pl.col('experiment_id'))

    for column in _gene_columns:
        if column in df.columns:
            expressions.append(pl.col(column))
        else:
            expressions.append(pl.lit(None).alias(column))

    extra_columns = [column for column in df.columns if column not in _gene_columns and column != 'experiment_id']
    if extra_columns:
        if len(extra_columns) == 1 and extra_columns[0] == 'other_info':
            expressions.append(pl.col('other_info'))
        else:
            struct_columns = [column for column in extra_columns if column != 'other_info']
            if 'other_info' in df.columns:
                struct_columns.append('other_info')
            expressions.append(pl.struct(struct_columns).struct.json_encode().alias('other_info'))
    elif 'other_info' in df.columns:
        expressions.append(pl.col('other_info'))
    else:
        expressions.append(pl.lit(None).alias('other_info'))

    return df.select(expressions)


def _to_metadata(metadata, ids, tables):
    metadata = list(metadata)
    meta_list = []
    to_remove_tables = []
    to_remove_metadata = []

    for table in tables:
        if isinstance(table, de_arrows):
            to_remove_tables.append(table)
            for item in table.id:
                index = metadata.index(table.experiment_metadata[item]) + table.id.index(item)
                meta_list.append(metadata[index])
                to_remove_metadata.append(metadata[index])

    for table in to_remove_tables:
        tables.remove(table)
    for item in to_remove_metadata:
        metadata.remove(item)

    for meta, table in zip(metadata, tables):
        if isinstance(table, de_arrow):
            meta_list.append(meta)
        else:
            m, extra_info = ExperimentMetadata().to_dict(meta, table)
            extra_info = json.loads(extra_info)
            for key, value in extra_info.items():
                m[key] = value
            meta_list.append(m)

    metadata_dict = {}
    for meta, item in zip(meta_list, ids):
        metadata_dict[item] = meta
    return metadata_dict




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
            df = _heal_genes(df, species)
        df = df.insert_column(0, pl.lit(experiment_id).alias('experiment_id'))
       
        return _clean_df(df), metadata_fields


   

    @classmethod
    def _from_arrow(cls, df, metadata, id, columns = None):
        """
        Convert a nanoarrow array back into a de_arrow object.
        
        Parameters
        ----------
        df : polars.DataFrame
            A polars DataFrame to convert
        metadata : dict
            A nested dictionary containing metadata {ID: {field: value}}
        
        Returns
        -------
        de_arrow
            A de_arrow object with the provided array and metadata
        """
        instance = cls.__new__(cls)
        instance._table = df
        instance.experiment_metadata = metadata
        instance.name = _get_experiment_attribute(metadata, 'experiment_name', id)
        instance.annotation_version = _get_experiment_attribute(metadata, 'annotation_version', id)
        instance.contrast = _get_experiment_attribute(metadata, 'contrast', id)
        instance.file = _get_experiment_attribute(metadata, 'file', id)
        instance.date = _get_experiment_attribute(metadata, 'date', id)
        instance.model = _get_experiment_attribute(metadata, 'model', id)
        instance.id = id
        return instance

        
    
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

    

    def __getattr__(self, name):
        table = object.__getattribute__(self, '_table')
        attr = getattr(table, name)
        if not callable(attr):
            return attr

        def wrapper(*args, **kwargs):
            result = attr(*args, **kwargs)
            if isinstance(result, pl.DataFrame):
                if result.columns == table.columns:
                    return self.__class__._from_arrow(result, self.experiment_metadata, self.id)
            
            return result

        return wrapper

    
    

    
    def get_significant_genes(self, log2fc=1, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        df = self._table.filter(
            (pl.col('log2fc').abs() >= log2fc) &
            (pl.col('pvalue') <= pvalue) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._from_arrow(df, self.experiment_metadata, self.id)

    def get_downregulated(self, log2fc=0, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        frame = self._table.filter(
            (pl.col('log2fc') <= log2fc) &
            (pl.col('pvalue') <= pvalue) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._from_arrow(df, self.experiment_metadata, self.id)

    def get_upregulated(self, log2fc=0, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        frame = self._table.filter(
            (pl.col('log2fc') >= log2fc) &
            (pl.col('pvalue') <= pvalue) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._from_arrow(frame, self.experiment_metadata, self.id)

    def set_id(self, id):
        if isinstance(id, int):
            frame = self._table.with_columns(pl.lit(id).alias('experiment_id'))
            metadata = {id: self.experiment_metadata[self.id]}
            return self._from_arrow(frame, metadata, id)
        raise DeQuackError(f'set_ids expected an integer, not {type(id)}')


    def get_gene(self, gene_symbol=None, ensembl_id=None, id=None):
        required_columns = {'gene_symbol', 'ensembl_id'}
        _check_columns(required_columns, self.columns)
        expression = pl.lit(True)
        if id is not None and 'experiment_id' in self.columns:
            expression = expression & (pl.col('experiment_id') == id)
        if gene_symbol is not None:
            expression = expression & (pl.col('gene_symbol') == gene_symbol)
        if ensembl_id is not None:
            expression = expression & (pl.col('ensembl_id') == ensembl_id)
        df = self._table.filter(expression)
        return self._from_arrow(df, self.experiment_metadata, self.id)
        

    def add_experiment(self, data, metadata = None, id=None):
        if isinstance(data, de_arrow):
            df, meta, new_ids = self.add_experiment_arrow(data, id)
        elif isinstance(data, de_arrows):
            df, meta, new_ids = self._add_experiment_arrows(data, id)
        else:
            df, meta, new_ids = self._add_experiment_data(data, metadata, id)
        return self._from_arrow(df, meta, new_ids)
    
    def add_experiment_arrow(self, data, id=None):
        df = self._table
        meta = self.experiment_metadata
        frame = data._table
        if id is not None:
            frame = frame.with_columns(pl.lit(id).alias('experiment_id'))
        df.extend(frame)
        id = id or data.id
        new_metadata = {id: data.experiment_metadata[data.id]}
        new_ids = [self.id, id]
        meta.update(new_metadata)
        return (df, meta, new_ids)

    def _add_experiment_arrows(self, data, id=None):
        df = self._table
        meta = self.experiment_metadata
        frame = data._table
        if id is not None:
            if isinstance(id, list):
                if len(id) != len(data.id):
                    raise DeQuackError('Number of provided ids does not match number of existing ids')
                id_dict = dict(zip(data.id, id))
                frame = frame.with_columns(pl.col('experiment_id').replace(id_dict).alias('experiment_id'))
                new_metadata = {n: value for n, value in zip(id, [data.experiment_metadata.values()])}
        new_metadata = data.experiment_metadata
        df = df.extend(frame)
        id = id or data.id
        new_ids = [self.id] + id
        meta.update(new_metadata)
        return (df, meta, new_ids)

    def _add_experiment_data(self, data, metadata, id):
        df = self._table
        meta = self.experiment_metadata
        frame = _to_polars_table(data)
        frame = _order_columns(frame)
        if id is None:
            raise DeQuackError('id must be provided when adding a non-de_arrow table')
        if not isinstance(id, int):
            raise DeQuackError('id must be an integer when adding a non-de_arrow table')
        if metadata is None:
            raise DeQuackError('metadata must be provided when adding a non-de_arrow table')
        if not isinstance(metadata, dict):
            raise DeQuackError('metadata must be a dictionary when adding a non-de_arrow table')
        frame = frame.insert_column(0, pl.lit(id).alias('experiment_id'))
        metadata_fields, extra_info=ExperimentMetadata().to_dict(metadata, data)
        for key, value in json.loads(extra_info).items():
            metadata_fields[key]=value
        new_metadata = {id: metadata_fields}
        frame = frame.with_columns(cs.numeric().fill_null(0), cs.string().fill_null(''), cs.struct().fill_null({}))
        df.extend(frame)
        new_ids = [self.id, id]
        meta.update(new_metadata)
        return (df, meta, new_ids)


    def insert(self, file, intialize_gene_table = False, species = None):
        required_columns = {'padj', 'pvalue', 'log2fc', 'gene_symbol', 'ensembl_id', 'logCPM', 'stat', 'other_info', 'experiment_id'}
        _check_columns(required_columns, self.columns)
        if not isinstance(file, str):
            raise TypeError(f'file must be a string, not {type(file)}')
        if not os.path.exists(os.path.abspath(file)):
            raise FileNotFoundError(f'File not found: {file}')
        file = os.path.abspath(file)
        df = self._table.select(pl.all().exclude(['experiment_id']))
        df = df.with_columns(
        pl.col("gene_symbol").replace("", None),
        pl.col("ensembl_id").replace("", None)
        )        
        metadata = self.experiment_metadata[self.id]
        metadata_fields, extra_info=ExperimentMetadata().to_dict(metadata, file)
        db = de_quackling(file).connect()
        if intialize_gene_table == True:
            if species is None:
                species = 'human'
                try:
                    db.intialize_gene_table(species)
                except DuplicateGeneTableError:
                    pass
        sample_data = db.conn.execute('SELECT * FROM (SELECT * FROM df LIMIT 100)').fetchall()
        data_signature = hashlib.md5(str(sample_data).encode()).hexdigest()
        result = db.conn.execute(
            "INSERT INTO experimental_data (model, date, file, experiment_name, contrast, annotation_version, normalization, other_info, data_signature) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?) RETURNING experiment_id",
            (metadata_fields.get('model'), metadata_fields.get('date'), metadata_fields.get('file'), metadata_fields.get('experiment_name'), metadata_fields.get('contrast'), metadata_fields.get('annotation_version'), metadata_fields.get('normalization'), extra_info, data_signature)
        ).fetchall()
        db.conn.execute(_core_queries['insert_de_arrow'], (result[0][0], species))
        

        
    def df(self):
        return self._frame.clone()


class de_arrows:
    def __new__(cls, *args, metadata = None, ids=None, keep_ids = False, heal_genes = False, species = None):
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
        if metadata is None:
            if len(none_de_arrow) == 0:
                metadata = []
            else:
                raise DeQuackError('Metadata must be provided for non-de_arrow tables')
        if len(args) < 1:
            raise DeQuackError('At least one table must be provided')
        if not isinstance(metadata, list):
            if len(none_de_arrow) == 1:
                metadata = [metadata]
            else:
                raise DeQuackError('Metadata must be provided as a list of dictionaries')
        if len(metadata) != len(none_de_arrow):
            raise DeQuackError('Number of metadata dictionaries does not match number of tables')
        if keep_ids == True:
            if ids is None:
                ids = [i + 1 for i in range(len(none_de_arrow))]
            if len(ids) != len(none_de_arrow):
                raise DeQuackError('Number of ids provided does not match number of tables')
            if isinstance(ids, int) and len(none_de_arrow) == 1:
                ids = [ids]
        else:
            if ids is None:
                if any(isinstance(arg, de_arrows) for arg in args):
                    raise DeQuackError('ids must be provided when combining de_arrows objects')
                ids = [i + 1 for i in range(len(args))]
            if len(ids) != len(args):
                raise DeQuackError('Number of ids provided does not match number of tables')
            
        _check_ids(ids)
        
        
        self._table, self.experiment_metadata, ids = self.__class__._from_tables(*args, metadata=metadata, ids=ids, heal_genes = heal_genes, species = species, keep_ids = keep_ids)
        self.id = ids
        self.name = _get_experiment_attribute(self.experiment_metadata, 'experiment_name', ids)
        self.annotation_version = _get_experiment_attribute(self.experiment_metadata, 'annotation_version', ids)
        self.contrast = _get_experiment_attribute(self.experiment_metadata, 'contrast', ids)
        self.file = _get_experiment_attribute(self.experiment_metadata, 'file', ids)
        self.date = _get_experiment_attribute(self.experiment_metadata, 'date', ids)
        self.model = _get_experiment_attribute(self.experiment_metadata, 'model', ids)

    
    def __getattr__(self, name):
        table = object.__getattribute__(self, '_table')
        attr = getattr(table, name)
        if not callable(attr):
            return attr

        def wrapper(*args, **kwargs):
            result = attr(*args, **kwargs)
            if isinstance(result, pl.DataFrame):
                if result.columns == table.columns:
                    return self.__class__._finalize_table(result, self.experiment_metadata, self.id)
            
            return result

        return wrapper
    
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
        args = list(args)
        for arg in final_args:
            if isinstance(arg, de_arrow) or isinstance(arg, de_arrows):
                for value in arg.experiment_metadata.values():
                    metadata.insert(final_args.index(arg), value)
                if keep_ids == True:
                    table.extend(arg._table)
                    final_ids.insert(final_args.index(arg), arg.id)
                    args.remove(arg)
        for arg, id in zip(final_args, final_ids):
            if isinstance(arg, de_arrow):
                arg = arg.with_columns(pl.lit(id).alias('experiment_id'))
                table.extend(arg._table)
            elif isinstance(arg, de_arrows):
                if not isinstance(id, list) and not isinstance(id, tuple) or len(id) != len(arg.id):
                    raise DeQuackError('ids must be a list or tuple for de_arrows objects and must be of same id length for de_arrows')
                id_lookup = pl.DataFrame({'old_id': arg.id, 'new_id': id}, schema = {'old_id': pl.Int32(), 'new_id': pl.Int32()})
                arg = arg.join(id_lookup, left_on = 'experiment_id', right_on = 'old_id', how = 'left')
                arg = arg.with_columns(pl.col('new_id').alias('experiment_id')).select(pl.all().exclude(['old_id', 'new_id']))
                table.extend(arg)
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
        table = table.with_columns(cs.numeric().fill_null(0), cs.string().fill_null(''), cs.struct().fill_null({}))
        return _clean_df(table), metadata_fields, flat

    @classmethod
    def _from_arrow(cls, table, metadata, ids, columns = None):
        instance = object.__new__(cls)
        instance._table = table
        instance.experiment_metadata = metadata
        instance.id = ids
        instance.name = _get_experiment_attribute(metadata, 'experiment_name', ids)
        instance.annotation_version = _get_experiment_attribute(metadata, 'annotation_version', ids)
        instance.contrast = _get_experiment_attribute(metadata, 'contrast', ids)
        instance.file = _get_experiment_attribute(metadata, 'file', ids)
        instance.date = _get_experiment_attribute(metadata, 'date', ids)
        instance.model = _get_experiment_attribute(metadata, 'model', ids)
        return instance

    def set_id(self, ids):
        if isinstance(ids, dict):
            return self._set_id_dict(ids)
        elif isinstance(ids, list):
            if len(ids) != len(self.id):
                raise DeQuackError('Number of provided ids does not match number of existing ids')
            return self._set_id_list(ids)
        raise DeQuackError(f'set_ids expected a dictionary or list, not {type(ids)}')
    
    def _set_id_dict(self, ids):
        df = self._table
        for old_id, new_id in ids.items():
            if old_id not in self.id:
                raise DeQuackError(f'{old_id} is not in existing ids. Existing ids are {self.id}')
        _check_ids(list(ids.values()) + self.id)
        df = df.with_columns(pl.col('experiment_id').replace(ids).alias('experiment_id'))
        new_ids = list(ids.values()) + [i for i in self.id if i not in ids.keys()]
        metadata = {new_id: self.experiment_metadata.get(old_id, {}) for old_id, new_id in ids.items()}
        for id in new_ids:
            if id not in metadata.keys():
                metadata[id] = self.experiment_metadata.get(id, {})

        return self._from_arrow(df, metadata, ids = new_ids)

    def _set_id_list(self, ids):
        if len(ids) != len(self.id):
            raise DeQuackError('Number of provided ids does not match number of existing ids')
        _check_ids(ids)
        mapping = dict(zip(self.id, ids))
        df = self._table.with_columns(pl.col('experiment_id').replace(mapping).alias('experiment_id'))
        metadata = {new_id: self.experiment_metadata.get(old_id, {}) for old_id, new_id in mapping.items()}
        return self._from_arrow(df, metadata, ids=ids)
    
    def get_experiment(self, id=None, name=None, model=None, annotation_version=None, normalization=None, date=None, contrast=None, file=None):

        selected_ids = set(self.id)
        if id is not None:
            if isinstance(id, list):
                selected_ids = {item for item in self.id if item in id}
            else:
                selected_ids = {item for item in self.id if item == id}


            def matches(meta):
                def text_match(field_value, expected):
                    if expected is None:
                        return True
                    if field_value is None:
                        return False
                    if isinstance(expected, list):
                        expected_values = [str(item).lower().strip() for item in expected]
                        return str(field_value).lower().strip() in expected_values
                    return str(expected).lower().strip() in str(field_value).lower().strip()

                return (
                    text_match(meta.get('experiment_name'), name) and
                    text_match(meta.get('model'), model) and
                    text_match(meta.get('annotation_version'), annotation_version) and
                    text_match(meta.get('normalization'), normalization) and
                    text_match(meta.get('contrast'), contrast) and
                    text_match(meta.get('file'), file) and
                    (date is None or str(meta.get('date')) == str(date))
                )

            selected_ids = {item for item in selected_ids if matches(self.experiment_metadata.get(item, {}))}

        if not selected_ids:
            raise DeQuackError('No experiments found with the provided metadata')

        frame = self._table.filter(pl.col('experiment_id').is_in(list(selected_ids))) if 'experiment_id' in self.columns else self._frame
        return self._finalize_table(frame)
    
    def _finalize_table(self, df):
        df_ids = df.select('experiment_id').unique().to_series().to_list() if 'experiment_id' in df.columns else []
        metadata = {}
        for id in df_ids:
            metadata[id] = self.experiment_metadata.get(id, {})
        if len(df_ids) < 2:
            return de_arrow._from_arrow(df, metadata, df_ids[0])
        return self._from_arrow(df, metadata, df_ids)
    
    def get_significant_genes(self, log2fc=1, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        df = self._table.filter(
            (pl.col('log2fc').abs() >= log2fc) &
            (pl.col('pvalue') <= pvalue) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._finalize_table(df)

    def get_downregulated(self, log2fc=0, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        df = self._table.filter(
            (pl.col('log2fc') <= log2fc) &
            (pl.col('pvalue') <= pvalue) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._finalize_table(df)

    def get_upregulated(self, log2fc=0, pvalue=0.05, logCPM=0):
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        df = self._table.filter(
            (pl.col('log2fc') >= log2fc) &
            (pl.col('pvalue') <= pvalue) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._finalize_table(df)
    



