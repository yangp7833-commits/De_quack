# DE_utilities
#
Python/R library powered by Python and DuckDB. Aims to help store differential analysis data using duckDB to improve reproductibility.
Try running this code to get started:

with DBManager('file.duckdb') as db:          
    db.initialize_gene_results('human')        
    db.insert_to_database('data.txt')        
    df=db.query(table='gene_results')        
    print(df)
