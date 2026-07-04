#!/home/codespace/.python/current/bin/python3
import pandas as pd
from de_quack import de_quackling, de_arrow, de_arrows
import duckdb
from duckdb import SQLExpression
import nanoarrow as na
df=pd.DataFrame([{'t':1, 'd':2}, {'t':4, 'd':5}])
metadata={'annotation_version':'ensembl_2', 'experiment_name':'test', 'contrast':'test', 't':'t'}
with de_quackling() as db:

    arrow = de_arrow('data1.txt', metadata = metadata)
    arrow.insert('SQL.duckdb')
    db.delete_experiment(id = 2)












