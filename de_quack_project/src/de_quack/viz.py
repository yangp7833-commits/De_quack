from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import re
import os
from .exceptions import ProcessingError
from .arrow import de_arrow, de_arrows, get_unique_conn
from .core import de_quackling
from .utilities import DE_ARROW_QUERIES


_de_queries = DE_ARROW_QUERIES

def volcano_plot(data, padj = 0.05, log2fc = 1, title=None, label_genes = 0, label_color = 'black', upregulated_color = 'red', downregulated_color = 'blue', insignificant_color = 'grey', label_size = 10, label_font = 'Arial', label_rotation = 0, label_fontweight = 'normal', file = None, label_name = 'gene_symbol'):
    de = de_quackling(get_unique_conn()).connect()
    if not isinstance(data, de_arrow) or not isinstance(data, de_arrows):
        de._preprocess(data)
        de_arrow_insertion_view = de._create_temp_view()
        arrow_table = de.conn.sql(_de_queries['insert_to_de_arrow'], params = {'id': 0})
    else:
        if data.columns - set(['gene', 'log2fc', 'padj', 'logCPM']) != set():
            raise ProcessingError("The input data must contain the following columns: 'gene', 'log2fc', and 'padj'")
        arrow_table = data

    
    upregulated = de.conn.sql('SELECT * FROM arrow_table WHERE log2fc > $log2fc AND padj < $padj AND padj != 0', params = {'padj': padj, 'log2fc': log2fc})
    insignificant = de.conn.sql('SELECT * FROM arrow_table WHERE abs(log2fc) < abs($log2fc) AND padj IS NOT NULL AND $padj != 0 OR padj > $padj AND padj != 0', params = {'padj': padj, 'log2fc': log2fc})
    downregulated = de.conn.sql('SELECT * FROM arrow_table WHERE log2fc < $log2fc AND padj < $padj AND padj != 0', params = {'padj': padj, 'log2fc': -1 * log2fc})

    log2fc_plot = [t[0] for t in upregulated.select('log2fc').fetchall()]
    padj_plot = [t[0] for t in upregulated.select('padj').fetchall()]
    plt.scatter(log2fc_plot, -np.log10(padj_plot), color=upregulated_color, alpha=0.5, label='Upregulated')
    log2fc_plot = [t[0] for t in downregulated.select('log2fc').fetchall()]
    padj_plot = [t[0] for t in downregulated.select('padj').fetchall()]
    plt.scatter(log2fc_plot, -np.log10(padj_plot), color=downregulated_color, alpha=0.5, label='Downregulated')
    log2fc_plot = [t[0] for t in insignificant.select('log2fc').fetchall()]
    padj_plot = [t[0] for t in insignificant.select('padj').fetchall()]
    plt.scatter(log2fc_plot, -np.log10(padj_plot), color=insignificant_color, alpha=0.5, label='Insignificant')
    plt.axhline(-np.log10(padj), color='black', linestyle='--', linewidth=1)
    plt.axvline(log2fc, color='black', linestyle='--', linewidth=1)

    plt.xlabel('Log2FC')
    plt.ylabel('Log10(Padj)')
    plt.title(title)
    plt.legend()

    if label_genes > 0:
        if label_name not in upregulated.columns or label_name not in downregulated.columns:
            raise ProcessingError(f"The label_name '{label_name}' is not present in the data. Please provide a valid column name for labeling genes.")
        top_regulated = upregulated.select(label_name, 'log2fc', 'padj').order('log2fc DESC').limit(label_genes).fetchall()
        top_regulated += downregulated.select(label_name, 'log2fc', 'padj').order('log2fc ASC').limit(label_genes).fetchall()
        for gene, log2fc, padj in top_regulated:
            plt.text(log2fc, -np.log10(padj), gene, fontsize=label_size, color=label_color, rotation=label_rotation, fontweight=label_fontweight, fontname=label_font)
       

    if file:
        if not os.path.exists(os.path.abspath(file)):
            open(os.path.abspath(os.path.abspath(file)), 'a').close()
        plt.savefig(file, bbox_inches='tight', dpi=300)
    

        



