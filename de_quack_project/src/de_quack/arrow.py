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
from typing import Sequence, TypeAlias
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

ExperimentId: TypeAlias = int 
ExperimentMetadataField: TypeAlias = str | int | float | bool | None | dict[str, object] | list[object]
ExperimentMetadataRecord: TypeAlias = dict[str, ExperimentMetadataField]
ExperimentMetadataMap: TypeAlias = dict[ExperimentId, ExperimentMetadataRecord]
MetadataTableSource: TypeAlias = pl.DataFrame | pl.LazyFrame | str | dict[str, object] | list[dict[str, object]]


def get_unique_conn() -> str:
    unique_id = uuid.uuid4().hex[:8]
    return f':memory:de_quack_{unique_id}'

def _clean_df(df: pl.DataFrame) -> pl.DataFrame:
    df = df.with_columns(pl.col('gene_symbol').fill_null(pl.lit('')).alias('gene_symbol'),
                        pl.col('ensembl_id').fill_null(pl.lit('')).alias('ensembl_id'),
                        pl.col('log2fc').fill_null(pl.lit(0.0)).alias('log2fc'),
                        pl.col('pvalue').fill_null(pl.lit(1.0)).alias('pvalue'),
                        pl.col('logCPM').fill_null(pl.lit(0.0)).alias('logCPM'),
                        pl.col('other_info').fill_null(pl.lit('{}')).alias('other_info'),
                        pl.col('padj').fill_null(pl.lit(1.0)).alias('padj'))
    return df


def _to_polars_table(table: object) -> pl.DataFrame:
    if isinstance(table, pl.DataFrame):
        return table
    if isinstance(table, pl.LazyFrame):
        return table.collect()

    try:
        import pandas as pd
    except ImportError:
        pd = None

    if pd is not None and isinstance(table, pd.DataFrame):
        return pl.from_pandas(table, nan_to_null = True)

    if isinstance(table, str):
        path = os.path.abspath(table)
        if not os.path.exists(path):
            raise FileNotFoundError(f'File not found: {path}')
        if path.lower().endswith('.parquet'):
            return pl.read_parquet(path)
        try:
            with open(path, 'r', encoding='utf-8') as handle:
                first_line = handle.readline()
        except Exception as e:
            raise FileNotFoundError(f'Error reading file: {path}') from e
        separator = '\t' if '\t' in first_line else (';' if ';' in first_line else ',')
        return pl.read_csv(path, separator=separator, null_values = ['', 'NA', 'NaN', 'nan'])

    if isinstance(table, dict):
        return pl.DataFrame(table)

    if isinstance(table, list):
        return pl.DataFrame(table)

    return pl.DataFrame(table)


def _check_columns(required_columns: set[str], columns: Sequence[str]) -> None:
    missing_columns = required_columns - set(columns)
    if missing_columns:
        raise DeQuackError(f'De_arrow object is missing {missing_columns} for this function')


def _check_ids(ids: Sequence[ExperimentId]) -> None:
    for item in ids:
        if ids.count(item) > 1:
            raise DeQuackError(f'id {item} occurs multiple times in object ids')
        if isinstance(item, str):
            raise TypeError('String found in id list')


def _get_experiment_attribute(
    metadata: ExperimentMetadataMap,
    field: str,
    ids: ExperimentId | Sequence[ExperimentId],
) -> dict[ExperimentId, ExperimentMetadataField | None]:
    if isinstance(ids, list):
        return {item: metadata.get(item, {}).get(field) for item in ids}
    return {ids: metadata.get(ids, {}).get(field)}


def _get_gene_table(species: str) -> object:
    folder = resources.files('de_quack') / 'gene_tables'
    return folder / f'{species}_genes.parquet'


def _heal_genes(df: pl.DataFrame, species: str) -> pl.DataFrame:
    exp_id = 'experiment_id' in df.columns
    df = df.lazy()
    genes = pl.scan_parquet(_get_gene_table(species)).select('symbol', 'ensembl_id')
    df = df.with_columns(pl.col('gene_symbol').fill_null(pl.lit('')).alias('gene_symbol'))
    df = df.join(genes, on='ensembl_id', how='left')
    df = df.with_columns(pl.col('ensembl_id').fill_null(pl.lit('')).alias('ensembl_id'))
    df = df.join(genes, left_on='gene_symbol', right_on='symbol', how='left', suffix='_on_symbol')
    if exp_id:
        df = df.select(
            pl.col('experiment_id'),
            pl.coalesce(['symbol', 'gene_symbol']).alias('gene_symbol'),
            pl.coalesce(['ensembl_id', 'ensembl_id_on_symbol']).alias('ensembl_id'),
            pl.all().exclude(['symbol', 'ensembl_id_on_symbol', 'gene_symbol', 'ensembl_id', 'experiment_id'])
        ).collect()
    else:
        df = df.select(
            pl.coalesce(['symbol', 'gene_symbol']).alias('gene_symbol'),
            pl.coalesce(['ensembl_id', 'ensembl_id_on_symbol']).alias('ensembl_id'),
            pl.all().exclude(['symbol', 'ensembl_id_on_symbol', 'gene_symbol', 'ensembl_id', 'experiment_id'])
        ).collect()
    return df


def _order_columns(df: object, columns: dict[str, str] | None = None) -> pl.DataFrame:
    if columns is None:
        columns = {}
    df = _to_polars_table(df)
    rename_map = {
        column: _GENE_ALIAS_TO_COLUMN[column.lower()]
        for column in df.columns
        if column.lower() in _GENE_ALIAS_TO_COLUMN and column != 'experiment_id'
    }
    for key, value in columns.items():
        if key in df.columns and key not in rename_map.keys():
            if value not in gene_columns:
                raise DeQuackError(f'Column {value} is not a valid gene column')
            rename_map[key] = value
        else:
            raise DeQuackError(f'Column {key} is not in the dataframe or is already mapped')
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
            if column == 'logCPM' and 'base_mean' in df.columns:
                expressions.append(pl.col('base_mean').add(1).log(base = 2).alias('logCPM'))
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






class DeArrow:

    def __init__(
        self,
        info: object,
        id: ExperimentId | None = None,
        metadata: ExperimentMetadataRecord | None = None,
        heal_genes: bool = False,
        species: str | None = None,
        columns: dict[str, str] | None = None,
        **fields: ExperimentMetadataField,
    ) -> None:
        if columns is None:
            columns = {}
        if metadata is None:
            metadata = fields
        if len(metadata) == 0:
            raise ProcessingError('No metadata provided in de_arrow object')
        
        if id is None:
            id = 1
        
        if not isinstance(metadata, dict):
            raise ProcessingError('Metadata must be a dictionary')
        self._table, self.experiment_metadata = self.__class__._to_de_arrow(info, metadata, columns, id, heal_genes = heal_genes, species = species)
        self.name = _get_experiment_attribute(self.experiment_metadata, 'experiment_name', id)
        self.annotation_version = _get_experiment_attribute(self.experiment_metadata, 'annotation_version', id)
        self.contrast = _get_experiment_attribute(self.experiment_metadata, 'contrast', id)
        self.file = _get_experiment_attribute(self.experiment_metadata, 'file', id)
        self.date = _get_experiment_attribute(self.experiment_metadata, 'date', id)
        self.model = _get_experiment_attribute(self.experiment_metadata, 'model', id)
        self.id = id
    
    @classmethod
    def _to_de_arrow(
        self,
        info: object,
        metadata: ExperimentMetadataRecord,
        columns: dict[str, str] | None,
        experiment_id: ExperimentId | None = None,
        heal_genes: bool = False,
        species: str | None = None,
    ) -> tuple[pl.DataFrame, ExperimentMetadataMap]:
        df = _to_polars_table(info)
        metadata_fields, extra_info=ExperimentMetadata().to_dict(metadata, info)
        for key, value in json.loads(extra_info).items():
            metadata_fields[key]=value
        metadata_fields={experiment_id: metadata_fields}
        df = _order_columns(df, columns)
        if heal_genes == True:
            if species is None:
                species = 'human'
            df = _heal_genes(df, species)
        df = df.insert_column(0, pl.lit(experiment_id).alias('experiment_id'))
       
        return _clean_df(df), metadata_fields


   

    @classmethod
    def _from_arrow(
        cls,
        df: pl.DataFrame,
        metadata: ExperimentMetadataMap,
        id: ExperimentId | None,
    ) -> "DeArrow":
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

        
    
    def __setattr__(self, name: str, value: object) -> None:
        if name in ('_table', 'experiment_metadata', 'name', 'annotation_version', 'contrast', 'file', 'date', 'id', 'model'):
            object.__setattr__(self, name, value)
        else:
            table = object.__getattribute__(self, '_table')
            setattr(table, name, value)
    
    def __getitem__(self, key: object) -> object:
        return self._table[key]
       

    def __len__(self) -> int:
        return len(self._table)
    
    def __repr__(self) -> str:
        return repr(self._table)
    
    def __str__(self) -> str:
        return(str(self._table))

    

    def __getattr__(self, name: str) -> object:
        table = object.__getattribute__(self, '_table')
        attr = getattr(table, name)
        if not callable(attr):
            return attr

        def wrapper(*args: object, **kwargs: object) -> object:
            result = attr(*args, **kwargs)
            if isinstance(result, pl.DataFrame):
                if result.columns == table.columns:
                    return self.__class__._finalize_table(self, result)
                else:
                    raise DeQuackError('DeArrows object columns are immutable.')
            return result

        return wrapper

    
    

    
    def get_significant_genes(self, log2fc: float = 1, pvalue: float = 0.05, logCPM: float = 0) -> "DeArrow":
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        df = self._table.filter(
            (pl.col('log2fc').abs() >= log2fc) &
            (pl.col('pvalue') <= pvalue) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._from_arrow(df, self.experiment_metadata, self.id)

    def get_downregulated(self, log2fc: float = 0, pvalue: float = 0.05, logCPM: float = 0) -> "DeArrow":
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        df = self._table.filter(
            (pl.col('log2fc') <= log2fc) &
            (pl.col('pvalue') <= pvalue) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._from_arrow(df, self.experiment_metadata, self.id)

    def get_upregulated(self, log2fc: float = 0, pvalue: float = 0.05, logCPM: float = 0) -> "DeArrow":
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        frame = self._table.filter(
            (pl.col('log2fc') >= log2fc) &
            (pl.col('pvalue') <= pvalue) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._from_arrow(frame, self.experiment_metadata, self.id)

    def set_id(self, id: ExperimentId) -> "DeArrow":
        if isinstance(id, int):
            frame = self._table.with_columns(pl.lit(id).alias('experiment_id'))
            metadata = {id: self.experiment_metadata[self.id]}
            return self._from_arrow(frame, metadata, id)
        raise DeQuackError(f'set_ids expected an integer, not {type(id)}')


    def get_gene(
        self,
        gene_symbol: str | None = None,
        ensembl_id: str | None = None,
        id: ExperimentId | None = None,
    ) -> "DeArrow":
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
        

    def add_experiment(
        self,
        data: object,
        metadata: ExperimentMetadataRecord | None = None,
        id: ExperimentId | None = None,
        **fields: ExperimentMetadataField,
    ) -> "DeArrows":

        ids = [self.id]
        if metadata is None and fields:
            metadata = fields 
        if isinstance(data, DeArrow):
            df, meta, new_ids = _add_experiment_arrow(self._table, self.experiment_metadata, ids, data, id)
        elif isinstance(data, DeArrows):
            df, meta, new_ids = _add_experiment_arrows(self._table, self.experiment_metadata, ids, data, id)
        else:
            df, meta, new_ids = _add_experiment_data(self._table, self.experiment_metadata, ids, data, metadata, id)
        return DeArrows._from_arrow(df, meta, new_ids)
    

    
    def insert(self, file: str, initialize_gene_table: bool = False, species: str | None = None) -> None:
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
        metadata_fields = self.experiment_metadata[self.id]
        extra_info = {}
        for key in metadata_fields.keys():
            if key not in ['model', 'date', 'file', 'experiment_name', 'contrast', 'annotation_version', 'normalization']:
                extra_info[key] = metadata_fields.pop(key)
        db = de_quackling(file).connect()
        if initialize_gene_table == True:
            species = species or 'human'
            try:
                db.initialize_gene_table(species)
            except DuplicateGeneTableError:
                raise DuplicateGeneTableError(f'Gene table for {species} already exists in database. Please use a different species or set initialize_gene_table to False')
        sample_data = df.select(pl.all()).limit(50)
        data_signature = hashlib.md5(str(sample_data).encode()).hexdigest()
        result = db.conn.execute(
            "INSERT INTO experimental_data (model, date, file, experiment_name, contrast, annotation_version, normalization, other_info, data_signature) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?) RETURNING experiment_id",
            (metadata_fields.get('model'), metadata_fields.get('date'), metadata_fields.get('file'), metadata_fields.get('experiment_name'), metadata_fields.get('contrast'), metadata_fields.get('annotation_version'), metadata_fields.get('normalization'), extra_info, data_signature)
        ).fetchall()
        db.conn.execute(_core_queries['insert_de_arrow'], (result[0][0], species))
        db.conn.commit()
        db.conn.close()
        

        
    def df(self) -> pl.DataFrame:
        return self._df.clone()



class DeArrows:
    def __new__(
        cls,
        *args: object,
        columns: dict[str, str] | None = None,
        metadata: list[ExperimentMetadataRecord] | ExperimentMetadataRecord | None = None,
        ids: list[ExperimentId] | None = None,
        keep_ids: bool = False,
        heal_genes: bool = False,
        species: str = 'human',
    ) -> object:
        if columns is None:
            columns = {}
        if metadata is None:
            return super().__new__(cls)
        if len(metadata) != len(args):
            if len(args) == 1 and isinstance(metadata, dict):
                return DeArrow(args[0], metadata=metadata)
        if len(args) == 1:
            return DeArrow(args[0], metadata=metadata[0])
        return super().__new__(cls)
    
    def __init__(
        self,
        *args: object,
        columns: dict[str, str] | None = None,
        metadata: list[ExperimentMetadataRecord] | ExperimentMetadataRecord | None = None,
        ids: list[ExperimentId] | None = None,
        heal_genes: bool = False,
        species: str = 'human',
        keep_ids: bool = False,
    ) -> None:
        if columns is None:
            columns = {}
        """Initialize a DeArrows object with multiple differential expression tables.
        
        Parameters
        ----------
        *args : str or file paths
            Variable length argument list of table references (file paths or identifiers)
        metadata : list of dict
            List of metadata dictionaries, one for each table.
        """

        
        none_de_arrow = [table for table in args if not isinstance(table, DeArrow) and not isinstance(table, DeArrows)]
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
                if any(isinstance(arg, DeArrows) for arg in args):
                    raise DeQuackError('ids must be provided when combining DeArrows objects')
                ids = [i + 1 for i in range(len(args))]
            if len(ids) != len(args):
                raise DeQuackError('Number of ids provided does not match number of tables')
            
        _check_ids(ids)
        
        
        self._table, self.experiment_metadata, ids = self.__class__._from_tables(*args, metadata=metadata, ids=ids, columns=columns, heal_genes = heal_genes, species = species, keep_ids = keep_ids)
        self.id = ids
        self.name = _get_experiment_attribute(self.experiment_metadata, 'experiment_name', ids)
        self.annotation_version = _get_experiment_attribute(self.experiment_metadata, 'annotation_version', ids)
        self.contrast = _get_experiment_attribute(self.experiment_metadata, 'contrast', ids)
        self.file = _get_experiment_attribute(self.experiment_metadata, 'file', ids)
        self.date = _get_experiment_attribute(self.experiment_metadata, 'date', ids)
        self.model = _get_experiment_attribute(self.experiment_metadata, 'model', ids)

    
    def __getattr__(self, name: str) -> object:
        table = object.__getattribute__(self, '_table')
        attr = getattr(table, name)
        if not callable(attr):
            return attr

        def wrapper(*args: object, **kwargs: object) -> object:
            result = attr(*args, **kwargs)
            if isinstance(result, pl.DataFrame):
                if result.columns == table.columns:
                    return self.__class__._finalize_table(self, result)
                else:
                    raise DeQuackError('DeArrows object columns are immutable.')
            return result

        return wrapper
    
    def __setattr__(self, name: str, value: object) -> None:
        if name in ('_table', 'experiment_metadata', 'name', 'annotation_version', 'contrast', 'file', 'date', 'id', 'model'):
            object.__setattr__(self, name, value)
        else:
            table = object.__getattribute__(self, '_table')
            setattr(table, name, value)
    
    def __getitem__(self, key: object) -> object:
        return self._table[key]
    
    def __len__(self) -> int:
        return len(self._table)
    
    def __repr__(self) -> str:
        return repr(self._table)
    
    def __str__(self) -> str:
        return(str(self._table))

    @classmethod
    def _from_tables(
        cls,
        *args: object,
        columns: dict[str, str] | None = {},
        metadata: list[ExperimentMetadataRecord] | None = None,
        ids: list[ExperimentId] | None = None,
        heal_genes: bool = False,
        species: str = 'human',
        keep_ids: bool = False,
    ) -> tuple[pl.DataFrame, ExperimentMetadataMap, list[ExperimentId]]:
        experiment_columns=['experiment_id', 'model', 'date', 'file', 'experiment_name', 'contrast', 'annotation_version', 'normalization', 'extra_info']
        table = pl.DataFrame(schema = {'experiment_id': pl.Int32(), 'gene_symbol': pl.String(), 'ensembl_id': pl.String(), 'log2fc': pl.Float64(), 'logCPM': pl.Float64(), 'pvalue': pl.Float64(), 'padj': pl.Float64(), 'stat': pl.Float64(), 'other_info': pl.String()})
        frames = []
        meta_by_id = {}
        flat_ids = []
        meta_iter, ids_iter = iter(metadata or [{}]), iter(ids or [])
        for arg in args:
            if isinstance(arg, DeArrow):
                new_id = arg.id if keep_ids == True else next(ids_iter)
                arg = arg.with_columns(pl.lit(new_id).alias('experiment_id'))
                frames.append(arg._table)
                meta_by_id[new_id] = arg.experiment_metadata[arg.id]
                flat_ids.append(new_id)
            elif isinstance(arg, DeArrows):
                new_ids = arg.id if keep_ids == True else next(ids_iter)
                if not isinstance(new_ids, list) and not isinstance(new_ids, tuple) or len(new_ids) != len(arg.id):
                    raise DeQuackError('ids must be a list or tuple for DeArrows objects and must be of same id length for DeArrows objects')
                if not keep_ids:
                    arg = arg.with_columns(pl.col('experiment_id').replace(dict(zip(arg.id, new_ids))).alias('experiment_id'))
                frames.append(arg._table)
                flat_ids.extend(new_ids)
                for old_id, new_id in zip(arg.id, new_ids):
                    meta_by_id[new_id] = arg.experiment_metadata[old_id]
            else:
                new_id = next(ids_iter)
                arg = _to_polars_table(arg)
                arg = _order_columns(arg, columns)
                arg = arg.insert_column(0, pl.lit(new_id).alias('experiment_id'))
                flat_ids.append(new_id)
                meta_by_id[new_id] = next(meta_iter)
                try:
                    frames.append(arg)
                except pl.exceptions.SchemaError as e:
                    raise DeQuackError(f'A table either has invalid data types or is missing required columns: {e}')
        table = pl.concat(frames, how = 'diagonal')
        if heal_genes == True:
            table = _heal_genes(table, species)
        _check_ids(flat_ids)
        return _clean_df(table), meta_by_id, flat_ids

    @classmethod
    def _from_arrow(
        cls,
        table: pl.DataFrame,
        metadata: ExperimentMetadataMap,
        ids: list[ExperimentId],
    ) -> "DeArrows":
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

    def set_id(self, ids: dict[ExperimentId, ExperimentId] | list[ExperimentId]) -> "DeArrows":
        if isinstance(ids, dict):
            return self._set_id_dict(ids)
        elif isinstance(ids, list):
            if len(ids) != len(self.id):
                raise DeQuackError('Number of provided ids does not match number of existing ids')
            return self._set_id_list(ids)
        raise DeQuackError(f'set_ids expected a dictionary or list, not {type(ids)}')
    
    def _set_id_dict(self, ids: dict[ExperimentId, ExperimentId]) -> "DeArrows":
        df = self._table
        for old_id, new_id in ids.items():
            if old_id not in self.id:
                raise DeQuackError(f'{old_id} is not in existing ids. Existing ids are {self.id}')
        _check_ids(list(ids.values()) + self.id)
        df = df.with_columns(pl.col('experiment_id').replace(ids).alias('experiment_id'))
        new_ids = self.id.copy()
        metadata = self.experiment_metadata.copy()
        for old_id, new_id in ids.items():
            index = self.id.index(old_id)
            new_ids[index] = new_id
            metadata[new_id] = self.experiment_metadata.get(old_id, {})

        return self._from_arrow(df, metadata, ids = new_ids)

    def _set_id_list(self, ids: list[ExperimentId]) -> "DeArrows":
        if len(ids) != len(self.id):
            raise DeQuackError('Number of provided ids does not match number of existing ids')
        _check_ids(ids)
        mapping = dict(zip(self.id, ids))
        df = self._table.with_columns(pl.col('experiment_id').replace(mapping).alias('experiment_id'))
        metadata = {new_id: self.experiment_metadata.get(old_id, {}) for old_id, new_id in mapping.items()}
        return self._from_arrow(df, metadata, ids=ids)
    
    def get_experiment(
        self,
        id: ExperimentId | list[ExperimentId] | None = None,
        name: str | list[str] | None = None,
        model: str | list[str] | None = None,
        annotation_version: str | list[str] | None = None,
        normalization: str | list[str] | None = None,
        date: str | None = None,
        contrast: str | list[str] | None = None,
        file: str | list[str] | None = None,
    ) -> "DeArrow | DeArrows":

        selected_ids = set(self.id)
        if id is not None:
            if isinstance(id, list):
                selected_ids = {item for item in self.id if item in id}
            else:
                selected_ids = {item for item in self.id if item == id}


            def matches(meta: ExperimentMetadataRecord) -> bool:
                def text_match(field_value: object, expected: object) -> bool:
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
    
    def _finalize_table(self, df: pl.DataFrame) -> "DeArrow | DeArrows":
        df_ids = df.select('experiment_id').unique().to_series().to_list() if 'experiment_id' in df.columns else []
        metadata = {}
        for id in df_ids:
            metadata[id] = self.experiment_metadata.get(id, {})
        if len(df_ids) == 0:
            return DeArrow._from_arrow(df, metadata, None)
        if len(df_ids) == 1:
            return DeArrow._from_arrow(df, metadata, df_ids[0])
        return self._from_arrow(df, metadata, df_ids)
    
    def get_significant_genes(self, log2fc: float = 1, padj: float = 0.05, logCPM: float = 0) -> "DeArrows":
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        df = self._table.filter(
            (pl.col('log2fc').abs() >= log2fc) &
            (pl.col('padj') <= padj) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._finalize_table(df)

    def get_downregulated(self, log2fc: float = -1, padj: float = 0.05, logCPM: float = 1) -> "DeArrows":
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        df = self._table.filter(
            (pl.col('log2fc') <= log2fc) &
            (pl.col('padj') <= padj) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._finalize_table(df)

    def get_upregulated(self, log2fc: float = 1, padj: float = 0.05, logCPM: float = 1) -> "DeArrows":
        required_columns = {'log2fc', 'pvalue', 'logCPM'}
        _check_columns(required_columns, self.columns)
        df = self._table.filter(
            (pl.col('log2fc') >= log2fc) &
            (pl.col('padj') <= padj) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._finalize_table(df)
    
    def get_gene(
        self,
        gene_symbol: str | None = None,
        ensembl_id: str | None = None,
        id: ExperimentId | None = None,
    ) -> "DeArrows":
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
        return self._finalize_table(df)
        
    
    def insert(self, file: str, initialize_gene_table: bool = False, species: str | None = None) -> None:
        required_columns = {'padj', 'pvalue', 'log2fc', 'gene_symbol', 'ensembl_id', 'logCPM', 'stat', 'other_info', 'experiment_id'}
        _check_columns(required_columns, self.columns)
        if not isinstance(file, str):
            raise TypeError(f'file must be a string, not {type(file)}')
        if not os.path.exists(os.path.abspath(file)):
            raise FileNotFoundError(f'File not found: {file}')
        file = os.path.abspath(file)
        df = self._table
        df = df.with_columns(
        pl.col("gene_symbol").replace("", None),
        pl.col("ensembl_id").replace("", None)
        )
        db = de_quackling(file).connect()
        old_ids = self.id
        new_ids = []
        for id in old_ids:        
            metadata_fields = self.experiment_metadata[id]
            extra_info = {}
            for key in metadata_fields.keys():
                if key not in ['model', 'date', 'file', 'experiment_name', 'contrast', 'annotation_version', 'normalization']:
                    extra_info[key] = metadata_fields.pop(key)
            sample_data = df.select(pl.all()).limit(50)
            data_signature = hashlib.md5(str(sample_data).encode()).hexdigest()
            result = db.conn.execute(
                "INSERT INTO experimental_data (model, date, file, experiment_name, contrast, annotation_version, normalization, other_info, data_signature) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?) RETURNING experiment_id",
                (metadata_fields.get('model'), metadata_fields.get('date'), metadata_fields.get('file'), metadata_fields.get('experiment_name'), metadata_fields.get('contrast'), metadata_fields.get('annotation_version'), metadata_fields.get('normalization'), extra_info, data_signature)
            ).fetchall()
            new_ids.append(result[0][0])
        id_map = pl.DataFrame({'old_id': old_ids, 'new_id': new_ids}, schema = {'old_id': pl.Int32(), 'new_id': pl.Int32()})
        if initialize_gene_table == True:
            species = species or 'human'
            try:
                db.initialize_gene_table(species)
            except DuplicateGeneTableError:
                raise DuplicateGeneTableError(f'Gene table for {species} already exists in database. Please use a different species or set initialize_gene_table to False')
        db.conn.execute(_core_queries['insert_de_arrows'], (species,))
        db.conn.commit()
        db.conn.close()
    
    def add_experiment(
        self,
        data: object,
        metadata: ExperimentMetadataRecord | None = None,
        id: ExperimentId | None = None,
        **fields: ExperimentMetadataField,
    ) -> "DeArrows":
        meta = dict(self.experiment_metadata)
        if metadata is not None and fields:
            metadata.update(fields)
        if isinstance(data, DeArrow):
            df, meta, new_ids = _add_experiment_arrow(self._table.clone(), meta, self.id, data, id)
        elif isinstance(data, DeArrows):
            df, meta, new_ids = _add_experiment_arrows(self._table.clone(), meta, self.id, data, id)
        else:
            df, meta, new_ids = _add_experiment_data(self._table.clone(), meta, self.id, data, metadata, id)
        return DeArrows._from_arrow(df, meta, new_ids)
    

def _add_experiment_arrow(
    df: pl.DataFrame,
    meta: ExperimentMetadataMap,
    ids: list[ExperimentId],
    data: DeArrow,
    id: ExperimentId | None = None,
) -> tuple[pl.DataFrame, ExperimentMetadataMap, list[ExperimentId]]:
    frame = data._table
    if id is not None:
        frame = frame.with_columns(pl.lit(id).alias('experiment_id'))
    df.extend(frame)
    id = id or data.id
    new_metadata = {id: data.experiment_metadata[data.id]}
    new_ids = ids + [id]
    meta.update(new_metadata)
    return (df, meta, new_ids)

def _add_experiment_arrows(
    df: pl.DataFrame,
    meta: ExperimentMetadataMap,
    ids: list[ExperimentId],
    data: DeArrows,
    new_id: ExperimentId | list[ExperimentId] | None = None,
) -> tuple[pl.DataFrame, ExperimentMetadataMap, list[ExperimentId]]:
    frame = data._table
    if new_id is not None:
        check_ids(ids)
        if isinstance(new_id, list):
            if len(new_id) != len(data.id):
                raise DeQuackError('Number of provided ids does not match number of existing ids')
            id_dict = dict(zip(data.id, new_id))
            frame = frame.with_columns(pl.col('experiment_id').replace(id_dict).alias('experiment_id'))
            new_metadata = {n: value for n, value in zip(new_id, data.experiment_metadata.values())}
    else:
        new_metadata = data.experiment_metadata
    df = df.extend(frame)
    new_id = new_id or data.id
    new_ids = ids + new_id
    meta.update(new_metadata)
    return (df, meta, new_ids)

def _add_experiment_data(
    df: pl.DataFrame,
    meta: ExperimentMetadataMap,
    ids: list[ExperimentId],
    data: object,
    metadata: ExperimentMetadataRecord,
    id: ExperimentId,
) -> tuple[pl.DataFrame, ExperimentMetadataMap, list[ExperimentId]]:
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
    frame = _clean_df(frame)
    df.extend(frame)
    new_ids = ids + [id]
    meta.update(new_metadata)
    return (df, meta, new_ids)
    

    



