#!/usr/bin/codespace/.python/current/bin/python3
from matplotlib import pyplot as plt
from adjustText import adjust_text
import pandas as pd
import numpy as np
from core import de_duckling
import os

def volcano_plot(df, padj=0.05, log2fc=1, plot_file=None, pvalue=None, significant_up='red', insignificant='grey', significant_down='blue', label=5):
    clean_headers = {
        'log2FoldChange': ['log2fc', 'log2fold', 'log2foldchange', 'logfc'], 
        'pvalue': ['p-value', 'pvalue'], 
        'padj': ['padj', 'false discovery rate', 'fdr'], 
        'gene_name': ['gene_name', 'gene'], 
        'logCPM': ['logcpm', 'basemean']
    }
    if isinstance(df, pd.DataFrame):
        pass
    elif isinstance(df, list) and isinstance(df[0], dict):
        df=pd.DataFrame(df)
    elif os.path.exists(os.path.abspath(df)):
        df=de_duckling().preprocess_df(df, id=None)
    else:
        print("Error: Input must be a DataFrame or a list of dictionaries or a csv/tsv file path")
        return
    rename_dict={original: standard for original in df.columns for standard, variants in clean_headers.items() if original.lower() in variants}

    if pvalue:
        column='pvalue'
        value=pvalue
    else:
        column='padj'
        value=padj
    df = df.rename(columns=rename_dict)
    df[column] = df[column].dropna()
    df[column] = df[column].replace(0, 1e-200)
    df['log2FoldChange'] = df['log2FoldChange'].dropna().replace(np.inf, 6).replace(-np.inf, -6)
    df['significant_up'] = (df[column] < value) & (df['log2FoldChange'] > log2fc)
    df['significant_down'] = (df[column] < value) & (df['log2FoldChange'] < -log2fc)

    significant_genes = df[df['significant_up'] | df['significant_down']]
    if label is not None and label > 0 and not significant_genes.empty:
        top_genes = significant_genes.sort_values('log2FoldChange', key=lambda x: x.abs(), ascending=False).head(label)
        

    plt.figure(figsize=(10, 6))
    plt.scatter(df['log2FoldChange'], -np.log10(df[column]), color=insignificant, alpha=0.5)
    plt.scatter(df[df['significant_up']]['log2FoldChange'], -np.log10(df[df['significant_up']][column]), color=significant_up, alpha=0.7)
    plt.scatter(df[df['significant_down']]['log2FoldChange'], -np.log10(df[df['significant_down']][column]), color=significant_down, alpha=0.7)
    plt.xlabel('Log2 Fold Change')
    plt.ylabel(f'-Log10 {column}')
    plt.title('Volcano Plot')
    plt.axhline(-np.log10(padj), color='black', linestyle='--')
    plt.axvline(log2fc, color='black', linestyle='--')
    plt.axvline(-log2fc, color='black', linestyle='--')

    if len(top_genes) > 0:
        texts = []
        for _, row in top_genes.iterrows():
            x = row['log2FoldChange']
            y = -np.log10(row[column])
            if row['gene_symbol'] is not None:
                texts.append(plt.text(x, y, str(row['gene_symbol']), fontsize=8, weight='bold'))
            else:
                texts.append(plt.text(x, y, str(row['ensembl_id']), fontsize=8, weight='bold'))
        adjust_text(texts, 
            # Repel from these coordinates
            x=x, y=y, 
            # Increase the repulsion area around points/text
            expand_points=(1.5, 1.5), 
            expand_text=(1.2, 1.2),
            # Add arrows to link labels to points
            arrowprops=dict(arrowstyle='-', color='black', lw=0.5),
            # Important for volcano plots: prevents labels from flying off the chart
            autoalign='xy',
            only_move={'points':'xy', 'text':'xy'})
    



    if plot_file:
        plt.savefig(plot_file)
    else:
        plt.show()
    
    return plt