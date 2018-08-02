import numpy as np, pandas as pd, itertools
from sklearn import mixture

## Open the dataframe, delete lncRNA without signal or too low signal in all epigenetic marker
X = pd.read_csv('../matrix/ten_cell_types_lncRNA.xls', sep='\t', index_col=0)

columns =[x for x in X.columns if x.find('single') == -1]

X = X[columns]

name = X['name']
del X['name']

if True:
    import numpy as np
    from sklearn.decomposition import PCA

    pca = PCA(n_components=6)
    pca.fit(X)

    explained = pca.explained_variance_ratio_
    print explained
    print np.sum(explained)
    print np.sum(explained[:6])
    new_X = pca.fit_transform(X)


lala
X['sum'] = X.sum(axis=1)



X = X[X['sum'] > 500]
del X['sum']


gmm = mixture.GaussianMixture(n_components=3, covariance_type='full')
gmm = gmm.fit(X)

labels = gmm.predict(X)

X['label'] = labels
X['name'] = name

X1 = X[X.label==0]
X2 = X[X.label==1]
X3 = X[X.label==2]

X1.to_csv('../matrix/X1_lncRNA.xls', sep='\t')
X2.to_csv('../matrix/X2_lncRNA.xls', sep='\t')
X3.to_csv('../matrix/X3_lncRNA.xls', sep='\t')



