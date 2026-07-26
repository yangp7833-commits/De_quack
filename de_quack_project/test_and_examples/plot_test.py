from de_quack import de_arrow, de_arrows, volcano_plot

def test_volcano_plot():
    arrow = de_arrow('data.txt', metadata = {'experiment_name': 'Test Experiment'}, heal_genes = True, species = 'human')
    plt = volcano_plot(arrow, title='Volcano Plot', label_genes = 5, label_color = 'black', label_type = 'ensembl_id', show=False)
    plt2 = volcano_plot('data1.txt', title='Volcano Plot 2', label_genes = 5, label_color = 'black', label_type = 'ensembl_id', show=False)