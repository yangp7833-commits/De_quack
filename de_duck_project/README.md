# de_duck

Python/R library for differential expression (DE) data management and visualization using DuckDB.

Stores DE analysis results with automatic gene annotation, enabling reproducible and efficient querying of experimental data.

## Features

- **Efficient Storage**: DuckDB-backed persistence for gene expression results
- **Automatic Annotation**: Maps gene symbols and Ensembl IDs with canonical references
- **Flexible Querying**: Filter results using simple keyword arguments (e.g., `padj__lt=0.05`)
- **Multi-language Support**: Native Python and R interfaces
- **Visualization**: Built-in volcano plot generation with publication-ready styling
- **Data Normalization**: Handles variable column names and formats automatically

## Installation

```bash
pip install de_duck
```

For Excel export support (optional):
```bash
pip install de_duck[excel]
```

## Quick Start

### Python

```python
from de_duck import de_duckling

with de_duckling('results.duckdb') as db:
    db.initialize_gene_table('human')
    db.insert_to_database('data.txt')
    results = db.query('gene_results', padj__lt=0.05)
    print(results)
    
    # Generate volcano plot
    db.volcano_plot(results, padj=0.05, log2fc=1, plot_file='volcano.png')
```

### R

```r
source('de_duck_project/R/wrapper.R')

duck <- de_duck()
duck$connect('results.duckdb')
duck$initialize_gene_table('human')
duck$insert_to_database('data.txt')
results <- duck$query('gene_results', padj__lt=0.05)
duck$volcano_plot(results, padj=0.05, log2fc=1, plot_file='volcano.png')
duck$close()
```

## Database Schema

### Tables

- **experimental_data**: Metadata about DE analysis experiments
- **gene_results**: Individual gene-level results with statistics
- **genes**: Reference gene annotations (human, etc.)

## License

MIT
