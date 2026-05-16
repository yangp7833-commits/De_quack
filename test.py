#!/home/codespace/.python/current/bin/python3 

from db_manager import DBManager
import re
import viz
with DBManager() as db:
  
  #db.insert_to_database('data.txt')
  df=db.query(table='gene_results')
  viz.volcano_plot('data.txt')
  
    
