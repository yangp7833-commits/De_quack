#!/home/codespace/.python/current/bin/python3
import pandas as pd
from de_quack import de_quackling, de_arrow, de_arrows, volcano_plot
import de_quack
import duckdb
from duckdb import SQLExpression
import nanoarrow as na
import time
import pyarrow as pa
def make_data(total):
    with open('new_data_csv', 'w') as f:
        i = 0
        f.write('gene,log2fc,padj,logCPM\n')
        while i < total:
            i+=1
            length = 11 - len(str(i))
            n = '0' * length + str(i)
            f.write(f'ENSG{n},2.5,0.01,5.0\n')
start_time = time.perf_counter()
arrow = de_arrow('data1.txt', experiment_name='new_data_csv').get_gene(ensembl_id = 'ENSG00000000001').set_id(3)
arrows = arrow.add_experiment(arrow)
print(arrows.experiment_metadata)
with de_quackling('SQL.duckdb') as db:
    db.delete_experiment(name = 'new_data_csv')
    arrows.insert('SQL.duckdb')
    df = db.get_gene()
    print(df)













