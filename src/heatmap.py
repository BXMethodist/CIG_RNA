import pandas as pd, os, numpy as np

X = pd.read_csv('../matrix/ten_cell_types_lncRNA_PCA.xls', sep='\t', index_col=0)
X['index'] = X.index
X['name'], X['cell_type'] = X['index'].str.split('_', 1).str
del X['index']

X1 = X.sort_values(by='PCA1', ascending=False).iloc[:1000, :]
X2 = X.sort_values(by='PCA2', ascending=False).iloc[:1000, :]

X1_lncRNAs = X1.name.unique()
X2_lncRNAs = X2.name.unique()

X1_lncData = {}
X2_lncData = {}

X1_candidates = X[X.name.isin(X1_lncRNAs)]
X2_candidates = X[X.name.isin(X2_lncRNAs)]

for x in X1_lncRNAs:
    cur_X = X1_candidates[X1_candidates.name == x]
    for i in cur_X.index:
        cell_type = cur_X.ix[i, 'cell_type']
        if cell_type not in X1_lncData.keys():
            X1_lncData[cell_type] = {}
        X1_lncData[cell_type][x] = cur_X.ix[i, 'PCA1']

for x in X2_lncRNAs:
    cur_X = X2_candidates[X2_candidates.name == x]
    for i in cur_X.index:
        cell_type = cur_X.ix[i, 'cell_type']
        if cell_type not in X2_lncData.keys():
            X2_lncData[cell_type] = {}
        X2_lncData[cell_type][x] = cur_X.ix[i, 'PCA2']

X1_df = pd.DataFrame(X1_lncData)
X1_df.index.name = 'gene_id'
X1_df.to_csv('PCA1_top.xls', sep='\t')

X2_df = pd.DataFrame(X2_lncData)
X2_df.index.name = 'gene_id'
X2_df.to_csv('PCA2_top.xls', sep='\t')


