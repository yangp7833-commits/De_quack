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
arrows = de_arrows('new_data_csv', arrow, metadata = {})
arrow = de_arrows(arrow, arrows, keep_ids = True)
print(arrow.id)
end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Time taken to create arrow: {elapsed_time:.4f} seconds")
print(arrows)











