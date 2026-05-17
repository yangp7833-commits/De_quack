from matplotlib import pyplot as plt
from adjustText import adjust_text
import pandas as pd
import numpy as np
import re
import os
from .exceptions import ProcessingError
import warnings
import logging
from contextlib import redirect_stdout



def _normalize_df_columns(df):
    """Normalize DataFrame column names to canonical forms for visualization.
    
    Accepts DataFrames with various column name conventions (e.g., 'log2FC', 'log2FoldChange', 
    'Log2fc', etc.) and maps them to standard names: log2fc, pvalue, padj, gene_symbol, 
    ensembl_id, logCPM.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input DataFrame with any column naming convention.
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with normalized column names and required columns created/filled.
    """
    insertion_columns = {
        'gene_symbol': ['symbol', 'gene_symbol', 'gene_name', 'gene'],
        'ensembl_id': ['ensembl_id', 'ensemblid', 'ensembl_gene_id', 'gene_id'],
        'log2fc': ['log2fc', 'log2foldchange', 'log2fold', 'logfc', 'log2fold_change', 'log_2_fold_change'],
        'logCPM': ['logcpm', 'basemean', 'basemean_log2'],
        'pvalue': ['pvalue', 'p-value', 'p_value', 'pval'],
        'padj': ['padj', 'fdr', 'false_discovery_rate', 'adjusted_pvalue', 'adj_pvalue'],
    }
    
    df = df.copy()
    
    # Build a case-insensitive flat mapping
    flat_map = {v.lower(): k for k, variants in insertion_columns.items() for v in variants}
    
    # Rename columns based on case-insensitive matching
    rename_dict = {col: flat_map[col.lower()] for col in df.columns if col.lower() in flat_map}
    df.rename(columns=rename_dict, inplace=True)
    
    # Ensure gene identifier columns exist and handle ambiguity
    if 'gene_symbol' not in df.columns and 'ensembl_id' not in df.columns:
        # If neither exists, create empty columns
        df['gene_symbol'] = None
        df['ensembl_id'] = None
    elif 'gene_symbol' in df.columns and 'ensembl_id' not in df.columns:
        # If only gene_symbol exists, create ensembl_id column
        df['ensembl_id'] = None
    elif 'ensembl_id' in df.columns and 'gene_symbol' not in df.columns:
        # If only ensembl_id exists, create gene_symbol column
        df['gene_symbol'] = None
    
    # Ensure all required columns exist
    required_cols = ['log2fc', 'pvalue', 'padj', 'logCPM', 'gene_symbol', 'ensembl_id']
    for col in required_cols:
        if col not in df.columns:
            df[col] = None
    
    return df


def volcano_plot(df, padj=0.05, log2fc=1, plot_file=None, pvalue=None, significant_up='red', 
                 insignificant='grey', significant_down='blue', label=5):
    """Generate a volcano plot from differential expression data.
    
    Parameters:
    -----------
    df : pd.DataFrame or list of dict
        DataFrame with gene expression data. Accepts flexible column naming conventions
        for log2fc, pvalue, padj, gene_symbol, and ensembl_id. Will automatically normalize
        common variants (e.g., 'log2FoldChange', 'Log2fc', 'logFC' all map to 'log2fc').
    padj : float
        Adjusted p-value threshold (default: 0.05).
    log2fc : float
        Log2 fold-change threshold (default: 1).
    plot_file : str, optional
        Path to save the plot. If None, displays interactively.
    pvalue : float, optional
        Use raw p-value instead of adjusted p-value if provided.
    significant_up : str
        Color for upregulated significant genes (default: 'red').
    insignificant : str
        Color for non-significant genes (default: 'grey').
    significant_down : str
        Color for downregulated significant genes (default: 'blue').
    label : int
        Number of top genes to label (default: 5).
    
    Returns:
    --------
    matplotlib.pyplot : pyplot module for further customization.
    """
    # Convert input to DataFrame
    if isinstance(df, pd.DataFrame):
        pass
    elif isinstance(df, list) and len(df) > 0 and isinstance(df[0], dict):
        df = pd.DataFrame(df)
    elif os.path.isfile(os.path.abspath(df)):
        df=pd.read_csv(df, sep=None, engine='Python')
    else:
        raise ProcessingError("Input must be a pandas DataFrame or a non-empty list of dictionaries or a file path.")
    
    # Normalize column names to standard format
    df = _normalize_df_columns(df)

    if pvalue:
        column = 'pvalue'
        value = pvalue
    else:
        column = 'padj'
        value = padj
    
    df[column] = df[column].fillna(1.0)
    df[column] = df[column].replace(0, 1e-200)
    df['log2fc'] = df['log2fc'].fillna(0).replace(np.inf, 6).replace(-np.inf, -6)
    df['significant_up'] = (df[column] < value) & (df['log2fc'] > log2fc)
    df['significant_down'] = (df[column] < value) & (df['log2fc'] < -log2fc)

    significant_genes = df[df['significant_up'] | df['significant_down']]
    top_genes = pd.DataFrame()
    if label is not None and label > 0 and not significant_genes.empty:
        top_genes = significant_genes.sort_values('log2fc', key=lambda x: x.abs(), ascending=False).head(label)

    plt.figure(figsize=(10, 6))
    plt.scatter(df['log2fc'], -np.log10(df[column]), color=insignificant, alpha=0.5)
    plt.scatter(df[df['significant_up']]['log2fc'], -np.log10(df[df['significant_up']][column]), color=significant_up, alpha=0.7)
    plt.scatter(df[df['significant_down']]['log2fc'], -np.log10(df[df['significant_down']][column]), color=significant_down, alpha=0.7)
    plt.xlabel('Log2 Fold Change')
    plt.ylabel(f'-Log10 {column}')
    plt.title('Volcano Plot')
    plt.axhline(-np.log10(padj), color='black', linestyle='--')
    plt.axvline(log2fc, color='black', linestyle='--')
    plt.axvline(-log2fc, color='black', linestyle='--')

    if not top_genes.empty:
        texts = []
        for _, row in top_genes.iterrows():
            x = row['log2fc']
            y = -np.log10(row[column])
            gene_label = str(row.get('gene_symbol', row.get('ensembl_id', '')))
            if gene_label and gene_label != 'None':
                texts.append(plt.text(x, y, gene_label, fontsize=8, weight='bold'))
        
        if texts:
            with warnings.catch_warnings():
                with redirect_stdout(open(os.devnull, 'w')):
                    warnings.filterwarnings("ignore")
                    adjust_text_logger = logging.getLogger('adjustText')
                    current_level = adjust_text_logger.getEffectiveLevel()
                    adjust_text_logger.setLevel(logging.ERROR)
                    try:
                        adjust_text(texts, 
                        expand_points=(1.5, 1.5), 
                        expand_text=(1.2, 1.2),
                        arrowprops=dict(arrowstyle='-', color='black', lw=0.5),
                        autoalign='xy',
                        only_move={'points': 'xy', 'text': 'xy'})
                    finally:
                        adjust_text_logger.setLevel(current_level)

    if plot_file:
        plt.savefig(plot_file)
    else:
        plt.show()
    
    return plt