#!/home/codespace/.python/current/bin/python3
import pandas as pd
from de_quack import de_quackling, de_arrow, de_arrows, volcano_plot
import duckdb
from duckdb import SQLExpression
import nanoarrow as na

df=pd.DataFrame([{'t':1, 'd':2}, {'t':4, 'd':5}])
metadata={'annotation_version':'ensembl_2', 'experiment_name':'test', 'contrast':'test', 't':'t'}
volcano_plot('data.txt', label_genes=2, file = 'plot.png')










