from de_quack import de_quackling
import os

class TestCore():
    def test_ingestion(self):
        with de_quackling() as db:
            db.initialize_gene_table('human')
            db.initialize_gene_table('mouse')
            db.ingest('data.txt', metadata = {'experiment_name': 'Test Experiment', 'description': 'This is a test experiment.'}, species = 'human', columns = {'log2FoldChange': 'log2fc'})
            db.ingest('data1.txt', metadata = {'experiment_name': 'Test Experiment 2', 'description': 'This is another test experiment.'}, species = 'mouse', columns = {'log2FoldChange': 'log2fc'})
            assert db.get_gene('AKT1').height > 0 and db.get_gene(ensembl_id = 'ENSMUSG00000000017').height > 0
    
    def test_queries(self):
        with de_quackling() as db:
            assert db.get_downregulated().height == 15
            assert db.get_upregulated().height == 22
            assert db.get_significant_genes().height == 37
            assert db.get_experiment(name = 'Test Experiment').height > 0 
    
    def test_remove(self):
        with de_quackling() as db:
            db.delete_experiment(name = 'Test Experiment')
            db.delete_experiment(name = 'Test Experiment 2')
            assert db.get_experiment(name = 'Test Experiment').height == 0
            assert db.get_experiment(name = 'Test Experiment 2').height == 0
            
            
    def test_cleanup(self):
        os.remove('SQL.duckdb')

    
    
