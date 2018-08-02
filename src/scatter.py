import numpy as np, pandas as pd, matplotlib.pyplot as plt

if True:
    X = pd.read_csv('../matrix/ten_cell_types_lncRNA.xls', sep='\t', index_col=0)
    columns = [x for x in X.columns if x.find('single') == -1 and x.find('signal') == -1]
    X = X[columns]
    name = X['name']
    del X['name']

    plt.scatter(X['h3k4me3_qn_total_width'], X['h3k27me3_qn_total_width'])
    plt.savefig('../plots/lncRNA_k4m3vk27m3_width.png')


if False:
    X = pd.read_csv('../matrix_gene/ten_cell_types_gene.xls', sep='\t', index_col=0)
    columns = [x for x in X.columns if x.find('single') == -1 and x.find('signal') == -1 and x!='gene_length' and x != 'RNA_exp']
    X = X[columns]
    name = X['name']
    del X['name']
    plt.scatter(X['h3k4me3_qn_total_width'], X['h3k27me3_qn_total_width'])
    plt.savefig('../plots/gene_k4m3vk27m3_width.png')