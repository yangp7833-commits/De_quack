#!/home/codespace/.python/current/bin/python3 

from core import de_quackling
import viz
with de_quackling() as db:
  
  db.insert_to_database('data.txt')

  
  viz.volcano_plot('data.txt', plot_file='volcano_plot.png')
  db.delete_experiments(experiment_id=3)
  db.delete_experiments(experiment_id=4)
  df=db.query(table='gene_results', save_path='data.csv')
  assert len(df)==0
  
    
