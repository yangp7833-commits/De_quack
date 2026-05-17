# de_quack

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
pip install de_quack
```

For Excel export support (optional):
```bash
pip install de_quack[excel]
```

## Quick Start

### Python

```python
from de_quack import de_quackling, volcano_plot

with de_quackling('results.duckdb') as db:
    db.initialize_gene_table('human')
    db.insert_to_database('data.txt', experiment_name='test', comparison_label='test')
    results = db.query('gene_results', padj__lt=0.05)
    print(results)
    
# Generate volcano plot
volcano_plot(results, padj=0.05, log2fc=1, plot_file='volcano.png')
```

### R

```r
# Find and source the wrapper directly from the pip installation directory
py_pkg_path <- dirname(reticulate::import("de_quack")$`__file__`)
source(file.path(py_pkg_path, "wrapper.R"))

# Initialize
duck <- de_quack(db_path = "results.duckdb")
duck$connect()
duck$initialize_gene_table('human')
duck$insert_to_database('data.txt', experiment_name = 'test', comparison_label = 'test')
results <- duck$query('gene_results', padj__lt=0.05)
duck$volcano_plot(results, padj=0.05, log2fc=1, plot_file='volcano.png')
duck$close()
```

## Database Schema

### Tables

- **experimental_data**: Metadata about DE analysis experiments
- **gene_results**: Individual gene-level results with statistics
- **genes**: Reference gene annotations (human, etc.)

### Advanced Querying & JSON Fallback
De_quack supports intuitive field modifiers (`__lt`, `__gt`, `__le`, `__ge`, `__ne`). If a metric isn't a native column in the database, DeDuck automatically queries it out of custom nested JSON metadata dynamically:

```python
# Regular column filtering + automatic nested JSON key metadata extraction
with de_quackling as db:
    results = db.query('gene_results', padj__lt=0.05, custom_biotype__ne='pseudogene')
```


## License

MIT
