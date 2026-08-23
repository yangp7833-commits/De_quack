
import re
import os
from .exceptions import ProcessingError
from .arrow import DeArrow, DeArrows, _to_polars_table, _order_columns
import polars as pl

def volcano_plot(df: object, padj: float = 0.05, log2fc: float = 1, title=None, show=False, label_genes = 0, label_type = 'ensembl_id', insignificant_color='grey', upregulated_color='red', downregulated_color='blue', file: str | None = None, **labeling_kwargs) -> 'matplotlib.figure.Figure':
    """ 
        Make a volcano plot from a DeArrow or DeArrows or polars or pandas DataFrame. The DataFrame must contain the columns 'log2fc' and 'padj'.
        The cutoffs for significance can be specified using the padj and log2fc parameters. The default values are 0.05 for padj and 1 for log2fc.
        None-DeArrow objects will first be converted to a polars dataframe using DeArrow methods. The plot will show upregulated genes in red, downregulated genes in blue, and insignificant genes in grey. 
        These can be customized using positional arguments.
        The user can specify the number of top regulated genes to label on the plot, as well as the type of label to use (e.g., 'ensembl_id', 'gene_name', etc.). 
        Additional labeling parameters can be passed as keyword arguments.
        Uses matplotlib for plotting. If matplotlib is not installed, an ImportError will be raised with instructions to install it.

    """
    defaults = {
        'color': 'black',
        'fontsize': 10,
        'rotation': 0,
        'fontweight': 'normal',
        'fontname': 'sans-serif',
        }
    wrong_keys = [key for key in labeling_kwargs.keys() if key not in defaults.keys()]
    if wrong_keys:
        raise ValueError(f"Invalid arguments provided: {', '.join(wrong_keys)}. Valid arguments are: {', '.join(defaults.keys())}.")
    for key in [k for k in defaults.keys() if k not in labeling_kwargs.keys()]:
        labeling_kwargs[key] = defaults[key]
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError('Matplotlib is required for volcano plotting. Please install it using "pip install matplotlib" or "pip install de_quack[plotting]"')
    if not isinstance(df, (DeArrow, DeArrows)):
        df_copy = _to_polars_table(df)
        df_copy = _order_columns(df_copy)
    else:
        df_copy = df.df()
    fig, ax = plt.subplots(figsize=(8, 6))
    df_copy = df_copy.filter(pl.col('log2fc').is_not_null() & pl.col('padj').is_not_null())
    df_copy = df_copy.filter(pl.col('padj') > 0)
    upregulated_df = df_copy.filter((pl.col('padj') < padj) & (pl.col('log2fc') > log2fc))
    downregulated_df = df_copy.filter((pl.col('padj') < padj) & (pl.col('log2fc') < -log2fc))
    insignificant_df = df_copy.filter((pl.col('padj') >= padj) | ((pl.col('log2fc') >= -log2fc) & (pl.col('log2fc') <= log2fc)))
    ax.scatter(upregulated_df['log2fc'], -upregulated_df['padj'].log10(), color=upregulated_color, label='Upregulated', alpha=0.7)
    ax.scatter(downregulated_df['log2fc'], -downregulated_df['padj'].log10(), color=downregulated_color, label='Downregulated', alpha=0.7)
    ax.scatter(insignificant_df['log2fc'], -insignificant_df['padj'].log10(), color=insignificant_color, label='Insignificant', alpha=0.7)

    ax.set_xlabel('Log2FC')
    ax.set_ylabel('-Log10(Padj)')
    ax.set_title(title)
    ax.legend()

    if label_genes > 0:
        if label_type in df_copy.columns:
            label_name = label_type
        else:
            raise ValueError(f"Label column '{label_type}' not found in the provided DataFrame columns.")

        top_regulated = upregulated_df.sort('log2fc', descending=True).select(pl.col(label_name), pl.col('padj'), pl.col('log2fc')).limit(label_genes).to_dicts()
        top_regulated += downregulated_df.sort('log2fc', descending=False).select(pl.col(label_name), pl.col('padj'), pl.col('log2fc')).limit(label_genes).to_dicts()
        import math
        for row in top_regulated:
            gene = row[label_name]
            gene_log2fc = row['log2fc']
            gene_padj = row['padj']
            ax.text(gene_log2fc + 0.01, -math.log10(gene_padj) + 0.1, gene, **labeling_kwargs)
       

    if file:
        if not os.path.exists(os.path.abspath(file)):
            os.makedirs(os.path.dirname(os.path.abspath(file)), exist_ok=True)
        fig.savefig(file, bbox_inches='tight', dpi=300)
    if show:
        fig.show()
    else:
        return fig
        
    

        



