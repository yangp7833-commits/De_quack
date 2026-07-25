from de_quack import de_arrow, de_arrows
import pytest
import os

class TestArrow():
    def test_arrow_creation(self):
        arrow_human = de_arrow('data.txt', metadata = {'experiment_name': 'Test Experiment'}, heal_genes = True, species = 'human')
        arrow_mouse = de_arrow('data1.txt', metadata = {'experiment_name': 'Test Experiment 2'})
        arrow_human_id = arrow_human.set_id(4)
        assert arrow_human._table.height > 0 and arrow_mouse._table.height > 0
        assert arrow_human.get_gene('AKT1').height > 0 and arrow_mouse.get_gene(ensembl_id = 'ENSMUSG00000000017').height > 0
        assert arrow_human_id.id == 4
    
    def test_arrows_creation(self):
        arrows = de_arrows('data.txt', 'data1.txt', heal_genes = True, species = 'human', metadata = [{'experiment_name1': 'Test Experiment'}, {'experiment_name': 'Test Experiment 2'}])
        arrows_with_ids = de_arrows(arrows, 'data.txt', metadata = {'experiment_name': 'Test Experiment 2'}, ids = [3], keep_ids = True)
        assert arrows._table.height > 0
        assert arrows.get_gene('AKT1').height > 0 and arrows.get_gene(ensembl_id = 'ENSMUSG00000000017').height > 0
        assert arrows_with_ids.id == [1, 2, 3]
    
    def test_arrow_queries(self):
        arrows = de_arrows('data.txt', 'data1.txt', metadata = [{'experiment_name': 'Test Experiment'}, {'experiment_name': 'Test Experiment 2'}], heal_genes = True)
        arrows_experiments = arrows.get_experiment(name = 'Test Experiment')
        assert arrows.get_downregulated().height == 15
        assert arrows.get_upregulated().height == 22
        assert arrows.get_significant_genes().height == 37
        assert len(arrows_experiments.experiment_metadata) == 2
    
    def test_arrow_add(self):
        arrow = de_arrow('data.txt', metadata = {'experiment_name': 'Test Experiment'}, heal_genes = True)
        arrows = arrow.add_experiment('data1.txt', metadata = {'experiment_name': 'Test Experiment 2'}, id = 3)
        assert arrows._table.height > 0
        assert arrows.get_gene('PTEN').height > 0 and arrows.get_gene(ensembl_id = 'ENSMUSG00000000017').height > 0
        assert arrows.id == [1, 3]
    
    def test_arrows_set_id(self):
        arrows = de_arrows('data.txt', 'data1.txt', metadata = [{'experiment_name': 'Test Experiment'}, {'experiment_name': 'Test Experiment 2'}])
        arrows1 = arrows.set_id([4, 5])
        arrows2 = arrows.set_id({2:7})
        assert arrows1.id == [4, 5]
        assert arrows2.id == [1, 7]
    

