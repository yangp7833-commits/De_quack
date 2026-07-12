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

#df = pd.read_csv('new_data_csv')
#df = df[df['log2fc'] > 1]
#print(df)
p = pa.Table.from_pydict({'g': [2, 3, 4]})


start_time = time.perf_counter()
arrow = de_arrows('new_data_csv', 'new_data_csv', ids =[1, 2], metadata = [{}] * 2)
print(arrow.experiment_metadata)
arrow = de_arrows._from_tables('new_data_csv', arrow, ids =[3], metadata = [{}], keep_ids = True)
end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Creating de_arrow took {elapsed_time:.2f} seconds.")
print(arrow)








