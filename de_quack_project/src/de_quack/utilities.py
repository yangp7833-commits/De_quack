from dataclasses import dataclass, field, fields, asdict
import difflib
import re
import os
from pathlib import Path
import datetime
from .exceptions import ProcessingError
import warnings
import logging
import json

_DATE_PATTERN = re.compile(r'(\d{4}-\d{2}-\d{2})')


# Configure logger with nice formatting
def _setup_logger():
    """Set up a logger with nice formatted output."""
    logger = logging.getLogger('de_quack')

    if not logger.handlers:
        logger.setLevel(logging.DEBUG)
        
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)

        formatter = logging.Formatter(
            '%(levelname)-8s | %(name)s | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    
    return logger


logger = _setup_logger()



def _warn(message: str, category=UserWarning, stacklevel=2):
    logger.warning(message)


def _info(message: str):
    logger.info(message)


def _debug(message: str):
    logger.debug(message)


def _error(message: str):
    """Log an error message."""
    logger.error(message)


DE_ARROW_QUERIES={


    'insert_to_de_arrow': 
                        '''
                        SELECT 
                                $id AS experiment_id,
                                COALESCE(ens.symbol, sym.symbol, e.gene_symbol) AS gene_symbol,
                                COALESCE(ens.ensembl_id, sym.ensembl_id, e.ensembl_id) AS ensembl_id,
                                COALESCE(TRY_CAST(e.log2fc AS DOUBLE), NULL) AS log2fc,
                                COALESCE(TRY_CAST(e.logCPM AS DOUBLE), log2(e.base_mean + 1), NULL) AS logCPM,
                                COALESCE(TRY_CAST(e.pvalue AS DOUBLE), NULL) AS pvalue,
                                COALESCE(TRY_CAST(e.padj AS DOUBLE), NULL) AS padj,
                                COALESCE(TRY_CAST(e.stat AS DOUBLE), NULL) AS stat,
                                e.other_info
                        FROM de_arrow_insertion_view as e
                            LEFT JOIN genes ens ON  (
                            e.ensembl_id IS NOT NULL AND ens.ensembl_id = UPPER(TRIM(e.ensembl_id))
                    )

                        LEFT JOIN genes sym ON (
                            e.gene_symbol IS NOT NULL AND sym.symbol = UPPER(TRIM(e.gene_symbol)) OR
                            e.gene_symbol IS NOT NULL AND sym.prev_symbol IS NOT NULL AND list_contains(sym.prev_symbol, e.gene_symbol)
                    )''',
        
    'map_mouse': '''SELECT 
                        *
                        FROM gene_table
                        LEFT JOIN csv_data ON gene_table.symbol = csv_data.column03 OR gene_table.ensembl_id = csv_data.column10
                        WITH csv_data AS (
                        SELECT 
                            csv_data.column03 as symbol,
                            CAST(regexp_replace(csv_data.column00, '^MGI:', '', 'i') AS INTEGER) AS id,
                            csv_data.column10 as ensembl_id
                        FROM
                        read_csv_auto('https://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt', strict_mode = false, delim = '\t', skip = 2))
                ''',
    'map_human': '''SELECT 
                csv_data.symbol as gene_symbol,
                csv_data.ensembl_gene_id as ensembl_id
                FROM read_csv_auto('https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt') AS csv_data
                '''
            
    }

DE_ARROWS_QUERIES={
    'insertion_from_table' : '''
                                WITH numbered_ids AS (
                                SELECT ids AS new_id, row_number() OVER () AS r_num FROM ids
                                    ),
                                numbered_data AS (
                                SELECT *, row_number() OVER () AS r_num FROM experimental_data)
                                    SELECT
                                        COALESCE(i.new_id, e.experiment_id) AS experiment_id,
                                        e.model,
                                        e.date,
                                        e.file,
                                        e.experiment_name,
                                        e.contrast,
                                        e.annotation_version,
                                        e.normalization,
                                        e.other_info AS extra_info,
                                        g.gene_symbol,
                                        g.ensembl_id,
                                        g.log2fc,
                                        g.logCPM,
                                        g.pvalue,
                                        g.padj,
                                        g.stat,
                                        g.other_info
                                FROM numbered_data e
                                LEFT JOIN numbered_ids i ON e.r_num = i.r_num
                                LEFT JOIN gene_results g
                                    ON e.experiment_id = g.experiment_id''',
    'set_ids' : '''             SELECT  
                                        COALESCE(i.new_id, g.experiment_id) AS experiment_id,      
                                        g.gene_symbol,
                                        g.ensembl_id,
                                        g.log2fc,
                                        g.logCPM,
                                        g.pvalue,
                                        g.padj,
                                        g.stat,
                                        g.other_info
                                FROM arrow_table g
                                LEFT JOIN ids_rel i ON i.old_id = g.experiment_id
                                ''',
                'get_experiment':   '''
                                    SELECT  
                                        g.experiment_id,      
                                        g.gene_symbol,
                                        g.ensembl_id,
                                        g.log2fc,
                                        g.logCPM,
                                        g.pvalue,
                                        g.padj,
                                        g.stat,
                                        g.other_info
                                FROM arrow_table g
                                WHERE
                                    g.experiment_id IN $ids
                                ''',
                'get_gene': '''
                            SELECT * FROM arrow_table
                            WHERE
                                ($id IS NULL OR experiment_id = $id) AND
                                ($gene_symbol IS NULL OR gene_symbol = $gene_symbol) AND
                                ($ensembl_id IS NULL OR ensembl_id = $ensembl_id)'''



    }

CORE_QUERIES = {
    'get_experiment': '''
                        SELECT 
                            e.experiment_id,
                            e.model,
                            e.date,
                            e.file,
                            e.experiment_name,
                            e.contrast,
                            e.annotation_version,
                            e.normalization,
                            e.other_info AS extra_info,
                        FROM experimental_data e
                        
                        WHERE 
                            ($1 IS NULL OR e.experiment_id = $1) AND
                            ($2 IS NULL OR e.date = $2) AND
                            ($3 IS NULL OR e.model LIKE '%' || $3 || '%') AND
                            ($4 IS NULL OR e.file LIKE '%' || $4 || '%') AND
                            ($5 IS NULL OR e.experiment_name LIKE '%' || $5 || '%') AND
                            ($6 IS NULL OR e.contrast LIKE '%' || $6 || '%') AND
                            ($7 IS NULL OR e.annotation_version LIKE '%' || $7 || '%') AND
                            ($8 IS NULL OR e.normalization LIKE '%' || $8 || '%')
                     ''',
    'get_significant': '''
                        SELECT
                            g.experiment_id,
                            g.gene_symbol,
                            g.ensembl_id,
                            g.log2fc,
                            g.logCPM,
                            g.pvalue,
                            g.padj,
                            g.stat,
                            g.other_info
                        FROM gene_results g
                        WHERE 
                            abs(g.log2fc) >= $log2fc AND
                            g.logCPM >= $logCPM AND
                            g.padj <= $padj AND
                            ($id IS NULL OR g.experiment_id = $id)
                            ''',
    'get_upregulated': '''
                        SELECT
                            g.experiment_id,
                            g.gene_symbol,
                            g.ensembl_id,
                            g.log2fc,
                            g.logCPM,
                            g.pvalue,
                            g.padj,
                            g.stat,
                            g.other_info
                        FROM gene_results g
                        WHERE 
                            g.log2fc >= $log2fc AND
                            g.logCPM >= $logCPM AND
                            g.padj <= $padj AND
                            ($id IS NULL OR g.experiment_id = $id)
                            ''',
    'get_downregulated': '''
                            SELECT
                            g.experiment_id,
                            g.gene_symbol,
                            g.ensembl_id,
                            g.log2fc,
                            g.logCPM,
                            g.pvalue,
                            g.padj,
                            g.stat,
                            g.other_info
                        FROM gene_results g
                        WHERE 
                            g.log2fc <= $log2fc AND
                            g.logCPM >= $logCPM AND
                            g.padj <= $padj AND
                            ($id IS NULL OR g.experiment_id = $id)''',
    
    'get_gene': '''
                    SELECT
                            g.experiment_id,
                            g.gene_symbol,
                            g.ensembl_id,
                            g.log2fc,
                            g.logCPM,
                            g.pvalue,
                            g.padj,
                            g.stat,
                            g.other_info
                        FROM gene_results g
                        WHERE 
                            ($id IS NULL OR g.experiment_id = $id) AND
                                ($gene_symbol IS NULL OR g.gene_symbol = $gene_symbol) AND
                            ($ensembl_id IS NULL OR g.ensembl_id = $ensembl_id)
                            ''',
    'find_delete_experiment': '''
                        SELECT experiment_id FROM experimental_data
                        WHERE
                            ($id IS NULL OR experiment_id = $id) AND
                            ($date IS NULL OR date = $date) AND
                            ($model IS NULL OR model LIKE '%' || $model || '%') AND
                            ($file IS NULL OR file LIKE '%' || $file || '%') AND
                            ($name IS NULL OR experiment_name LIKE '%' || $name || '%') AND
                            ($contrast IS NULL OR contrast LIKE '%' || $contrast || '%') AND
                            ($annotation_version IS NULL OR annotation_version LIKE '%' || $annotation_version || '%') AND
                            ($normalization IS NULL OR normalization LIKE '%' || $normalization || '%')
                        ''',
    'delete_experiment': ''' 
                        DELETE FROM experimental_data
                        WHERE 
                            experiment_id = ANY($ids)
                        ''',
    'delete_gene_results': '''
                        DELETE FROM gene_results
                        WHERE 
                            experiment_id = ANY($ids)
                        '''
    }

gene_mapping = {'mouse_genes': '''
            INSERT INTO genes (symbol, id, ensembl_id, alias_symbol, prev_symbol, species)
            SELECT 
                UPPER(TRIM(csv_data.column03)) as symbol,
                CAST(regexp_replace(csv_data.column00, '^MGI:', '', 'i') AS INTEGER) AS id, 
                csv_data.column10 as ensembl_id,
                NULL AS alias_symbol, 
                NULL AS prev_symbol,
                'mouse' as species 
            FROM read_csv_auto('https://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt', strict_mode = false, delim = '\t', skip = 2) AS csv_data
            ''',
            'human_genes': '''
            INSERT INTO genes (symbol, id, ensembl_id, alias_symbol, prev_symbol, species)
            SELECT 
                UPPER(TRIM(csv_data.symbol)) as symbol, 
                CAST(regexp_replace(csv_data.hgnc_id, '^hgnc:', '', 'i') AS INTEGER) AS id, 
                csv_data.ensembl_gene_id AS ensembl_id, 
                string_split(UPPER(COALESCE(csv_data.alias_symbol, '')), '|') as alias_symbol, 
                string_split(UPPER(COALESCE(csv_data.prev_symbol, '')), '|') as prev_symbol, 
                'human' as species 
            FROM read_csv_auto('https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt') AS csv_data
            '''
            }

gene_columns = {
        'gene_symbol': [
            'gene_symbol', 'symbol', 'hgnc_symbol', 'genesymbol',
        ],
        'gene_name':['gene', 'gene_name', 'genename', 
            'name'],
        'ensembl_id': [
            'ensembl_id', 'ensembl', 'gene_id', 'geneid', 'ensembl_gene_id', 
            'target_id', 'feature_id'
        ],
        'pvalue': [
        'pvalue', 'p_value', 'p.value', 'p-value', 'pval', 'p.val', 'p_val'
        ],
        'padj': [
        'padj', 'p.adj', 'p_adj', 'p.adjusted', 'p_adjusted', 
        'fdr', 'qval', 'qvalue', 'q-value', 'adj.p.val', 'adj.p.value'
        ],
        'stat': [
            'stat', 'f', 'lr', 't'
        ],
        'log2fc': [
            'log2foldchange', 'log2fc', 'log2_fc', 'log2.fc', 'logfc', 'log_fc', 'log.fc'
            ],
        'base_mean': [
        'basemean', 'base_mean',  'aveexpr', 'tpm', 'fpkm'
        ],
        'logCPM':['logcpm', 'log_cpm', 'log.cpm']
        }


def _find_date_and_file(info):
    if isinstance(info, str):
        info=os.path.abspath(info)
    else:
        return datetime.datetime.now().strftime('%Y-%m-%d'), 'in-memory object'

    if os.path.isfile(info):
        match = _DATE_PATTERN.search(info) # first tries to find the date using regex on the file name
        file_path=os.path.abspath(info)
        if match:
            return datetime.datetime.strptime(match.group(1), '%Y-%m-%d'), file_path
        else: # if the search doesn't work, then uses the path library to find the time, or else uses none
            path=Path(file_path)
            if path.stat().st_mtime:
                return datetime.datetime.fromtimestamp(path.stat().st_mtime).strftime('%Y-%m-%d'), file_path
            else:
                return datetime.datetime.now().strftime('%Y-%m-%d'), file_path


@dataclass
class ExperimentMetadata:
    experiment_name: str = None
    date: str = None
    contrast: str = None
    file: str = None
    normalization: str = None
    model: str = None
    annotation_version: str = None
    
    additional_info: dict = field(default_factory=dict)



    def to_dict(self, metadata: dict, info):
        """Convert the dataclass to a dictionary, including additional_info."""
        field_names = [f.name for f in fields(self)]
        field_name_set = set(field_names)
        for key, value in metadata.items():
            cleaned_key = key.lower().strip().replace(' ', '_')
            if cleaned_key in field_name_set:
                setattr(self, cleaned_key, value)
            else:
                matches = difflib.get_close_matches(cleaned_key, field_names, n=1, cutoff=0.8)
                if matches:
                    setattr(self, matches[0], value)
                else:
                    self.additional_info[key] = value

        core_dict = {
            'experiment_name': self.experiment_name,
            'date': self.date,
            'contrast': self.contrast,
            'file': self.file,
            'normalization': self.normalization,
            'model': self.model,
            'annotation_version': self.annotation_version,
            'additional_info': dict(self.additional_info),
        }
        if core_dict['date'] is None:
            self.date, _ = _find_date_and_file(info)
            core_dict['date'] = self.date
        else:
            try:
                core_dict['date'] = datetime.datetime.strptime(core_dict['date'], '%Y-%m-%d').strftime('%Y-%m-%d')
            except ValueError:
                _warn(f"Provided date '{core_dict['date']}' is not in the correct format (YYYY-MM-DD). Defaulting to current date.")
                self.date, _ = _find_date_and_file(info)
                core_dict['date'] = self.date
        if core_dict['file'] is None:
            _, file_path = _find_date_and_file(info)
            core_dict['file'] = file_path
        if core_dict['experiment_name'] is None and 'file' in self.additional_info:
            self.experiment_name = core_dict['file']
            core_dict['experiment_name'] = self.experiment_name 
        for field_name in field_names:
            if core_dict[field_name] is None:
                _warn(f'{field_name} is not provided in the metadata. {field_name} has been defaulted to unknown')
                core_dict[field_name] = 'unknown'

        other_info = core_dict.pop("additional_info")

        return core_dict, json.dumps(other_info)


       

