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
        
        'get_important_genes' : "SELECT * FROM arrow_table WHERE abs(log2fc) > $log2fc AND logCPM > $logCPM AND pvalue <= $pvalue",
        'get_downregulated' : "SELECT * FROM arrow_table WHERE log2fc < $log2fc AND logCPM > $logCPM AND pvalue <= $pvalue",
        'get_upregulated' : "SELECT * FROM arrow_table WHERE log2fc > $log2fc AND logCPM > $logCPM AND pvalue <= $pvalue",
        'set_arrow_id' : "SELECT * REPLACE ($id AS experiment_id) FROM arrow_table",
        'get_gene': 'SELECT * FROM arrow_table WHERE ($gene_symbol IS NULL OR gene_symbol = $gene_symbol) AND ($ensembl_id IS NULL OR gene_symbol = $gene_symbol)'

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
                            g.experiment_id,
                            g.gene_symbol,
                            g.ensembl_id,
                            g.log2fc,
                            g.logCPM,
                            g.pvalue,
                            g.padj,
                            g.stat,
                            g.other_info
                        FROM experimental_data e
                        LEFT JOIN gene_results g ON
                            e.experiment_id = g.experiment_id
                        
                        WHERE 
                            ($experiment_id IS NULL OR e.experiment_id = $experiment_id) AND
                            ($date IS NULL OR e.date = $date) AND
                            ($model IS NULL OR e.model LIKE '%' || $model || '%') AND
                            ($file IS NULL OR e.file LIKE '%' || $file || '%') AND
                            ($experiment_name IS NULL OR e.experiment_name LIKE '%' || $experiment_name || '%') AND
                            ($contrast IS NULL OR e.contrast LIKE '%' || $contrast || '%') AND
                            ($annotation_version IS NULL OR e.annotation_version LIKE '%' || $annotation_version || '%') AND
                            ($normalization IS NULL OR e.normalization LIKE '%' || $normalization || '%')
                     ''',
    'get_significant': '''
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
                            g.experiment_id,
                            g.gene_symbol,
                            g.ensembl_id,
                            g.log2fc,
                            g.logCPM,
                            g.pvalue,
                            g.padj,
                            g.stat,
                            g.other_info
                        FROM experimental_data e
                        LEFT JOIN gene_results g ON
                            e.experiment_id = g.experiment_id
                        WHERE 
                            abs(g.log2fc) >= $log2fc AND
                            g.logCPM >= $logCPM AND
                            g.padj <= $padj
                            ''',
    'get_upregulated': '''
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
                            g.experiment_id,
                            g.gene_symbol,
                            g.ensembl_id,
                            g.log2fc,
                            g.logCPM,
                            g.pvalue,
                            g.padj,
                            g.stat,
                            g.other_info
                        FROM experimental_data e
                        LEFT JOIN gene_results g ON
                            e.experiment_id = g.experiment_id
                        WHERE 
                            g.log2fc >= $log2fc AND
                            g.logCPM >= $logCPM AND
                            g.padj <= $padj''',
    'get_downregulated': '''
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
                            g.experiment_id,
                            g.gene_symbol,
                            g.ensembl_id,
                            g.log2fc,
                            g.logCPM,
                            g.pvalue,
                            g.padj,
                            g.stat,
                            g.other_info
                        FROM experimental_data e
                        LEFT JOIN gene_results g ON
                            e.experiment_id = g.experiment_id
                        WHERE 
                            g.log2fc <= $log2fc AND
                            g.logCPM >= $logCPM AND
                            g.padj <= $padj''',
    
    'get_gene': '''
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
                            g.experiment_id,
                            g.gene_symbol,
                            g.ensembl_id,
                            g.log2fc,
                            g.logCPM,
                            g.pvalue,
                            g.padj,
                            g.stat,
                            g.other_info
                        FROM experimental_data e
                        LEFT JOIN gene_results g ON
                            e.experiment_id = g.experiment_id
                        WHERE 
                            ($id IS NULL OR g.experiment_id = $id) AND
                            ($gene_symbol IS NULL OR g.gene_symbol = $id) AND
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
        date_pattern = re.compile(r'(\d{4}-\d{2}-\d{2})') 
        match = date_pattern.search(info) # first tries to find the date using regex on the file name
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
        for key, value in metadata.items():
            cleaned_key = key.lower().strip().replace(' ', '_')
            if cleaned_key in [f.name for f in fields(self)]:
                setattr(self, cleaned_key, value)
            elif difflib.get_close_matches(cleaned_key, [f.name for f in fields(self)], n=1, cutoff=0.8):
                closest_match = difflib.get_close_matches(cleaned_key, [f.name for f in fields(self)], n=1, cutoff=0.8)[0]
                setattr(self, closest_match, value)
            else:
                self.additional_info[key] = value
            
        core_dict = asdict(self)
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
        for field in [f.name for f in fields(self)]:
            if core_dict[field] is None:
                _warn(f'{field} is not provided in the metadata. {field} has been defaulted to unknown')
                core_dict[field] = 'unknown'
        
        if core_dict['additional_info'] is None:
            core_dict['additional_info'] = {}
        other_info = core_dict.pop("additional_info")
            
        return core_dict, json.dumps(other_info)


       

