import numpy as np, pandas as pd, itertools, matplotlib.pyplot as plt
from sklearn import mixture
from sklearn.decomposition import PCA


## Open the dataframe, delete lncRNA without signal or too low signal in all epigenetic marker
if False:
    from sklearn.decomposition import PCA

    X = pd.read_csv('../matrix/ten_cell_types_lncRNA.xls', sep='\t', index_col=0)

    columns = [x for x in X.columns if x.find('single') == -1 and x.find('signal') == -1]

    X = X[columns]

    name = X['name']
    del X['name']
    pca = PCA(n_components=6)
    pca.fit(X)

    new_X = pca.fit_transform(X)
    df = pd.DataFrame(new_X, index=X.index, columns=['PCA'+str(i) for i in range(1,7)])
    df.to_csv('../matrix/ten_cell_types_lncRNA_PCA.xls', sep='\t')

    plt.scatter(new_X[:, 0], new_X[:, 1])
    plt.savefig('../plots/ten_celltypes_PCA.png')

    cell_types = ['CD34', 'GM12878', 'H1-hESC', 'HMEC', 'HSMM', 'MRG', 'MSC', 'neural', 'NHLF', 'HUVEC']

    for cell_type in cell_types:
        X = pd.read_csv('../matrix/'+cell_type+'_parameters_real_table.xls', sep='\t', index_col=0)

        columns = [x for x in X.columns if x.find('single') == -1]

        X = X[columns]

        name = X['name']
        del X['name']

        new_X = pca.fit_transform(X)

        plt.scatter(new_X[:, 0], new_X[:, 1])
        plt.savefig('../plots/'+cell_type+'_PCA.png')


## get the components for PCA
if True:
    from sklearn.decomposition import PCA

    X = pd.read_csv('../matrix/ten_cell_types_lncRNA.xls', sep='\t', index_col=0)
    # X = pd.read_csv('../matrix_gene/ten_cell_types_gene.xls', sep='\t', index_col=0)

    columns = [x for x in X.columns if x.find('single') == -1 and x.find('signal') == -1 and x!='gene_length' and x != 'RNA_exp']
    print columns

    X = X[columns]

    name = X['name']
    del X['name']
    pca = PCA(n_components=6)
    pca.fit(X)

    print pca.explained_variance_ratio_

    df = pd.DataFrame(pca.components_, columns=list(X.columns))

    df.to_csv('../PCA_contributions.xls', sep='\t')
    # df.to_csv('../PCA_contributions_gene.xls', sep='\t')