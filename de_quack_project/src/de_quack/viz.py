
import re
import os
from .exceptions import ProcessingError
from .arrow import DeArrow, DeArrows, _to_polars_table, _order_columns
from .core import de_quackling
import polars as pl
from .utilities import DE_ARROW_QUERIES


_de_queries = DE_ARROW_QUERIES

def volcano_plot(df, padj = 0.05, log2fc = 1, title=None, label_genes = 0, label_color = 'black', upregulated_color = 'red', downregulated_color = 'blue', insignificant_color = 'grey', label_size = 10, label_font = 'Courier New', label_rotation = 0, label_fontweight = 'bold', file = None, label_type = 'ensembl_id', show=False):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError('Matplotlib is required for volcano plotting. Please install is using "pip install matplotlib" or "pip install de_quack[plotting]"')
    if not isinstance(df, (de_arrow, de_arrows)):
        df = _to_polars_table(df)
        df = _order_columns(df)
    df = df.filter(pl.col('log2fc').is_not_null() & pl.col('padj').is_not_null())
    df = df.filter(pl.col('padj') > 0)
    upregulated_df = df.filter((pl.col('padj') < padj) & (pl.col('log2fc') > log2fc))
    downregulated_df = df.filter((pl.col('padj') < padj) & (pl.col('log2fc') < -log2fc))
    insignificant_df = df.filter((pl.col('padj') >= padj) | ((pl.col('log2fc') >= -log2fc) & (pl.col('log2fc') <= log2fc)))
    plt.scatter(upregulated_df['log2fc'], -upregulated_df['padj'].log10(), color=upregulated_color, label='Upregulated', alpha=0.7)
    plt.scatter(downregulated_df['log2fc'], -downregulated_df['padj'].log10(), color=downregulated_color, label='Downregulated', alpha=0.7)
    plt.scatter(insignificant_df['log2fc'], -insignificant_df['padj'].log10(), color=insignificant_color, label='Insignificant', alpha=0.7)

    plt.xlabel('Log2FC')
    plt.ylabel('-Log10(Padj)')
    plt.title(title)
    plt.legend()

    if label_genes > 0:
        
        if label_type in upregulated_df.columns:
            label_name = label_type
        else:
            raise ValueError(f"Label column '{label_type}' not found in DataFrame columns.")

        top_regulated = upregulated_df.select(pl.col(label_name), pl.col('padj'), pl.col('log2fc').sort(descending=True)).limit(label_genes).to_dicts()
        top_regulated += downregulated_df.select(pl.col(label_name), pl.col('padj'), pl.col('log2fc').sort(descending=False)).limit(label_genes).to_dicts()
        for row in top_regulated:
            gene = row[label_name]
            log2fc = row['log2fc']
            padj = row['padj']
            plt.text(log2fc + 0.01, -pl.col('padj').log(base=10) + 0.1, gene, fontsize=label_size, color=label_color, rotation=label_rotation, fontweight=label_fontweight, fontname=label_font)
       

    if file:
        if not os.path.exists(os.path.abspath(file)):
            open(os.path.abspath(os.path.abspath(file)), 'a').close()
        plt.savefig(file, bbox_inches='tight', dpi=300)
    if show == True:
        plt.show()
    else:
        return plt
        
    

        



