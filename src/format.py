import pandas as pd, numpy as np, os

files = [x for x in os.listdir('../results/') if x.find('gene')!=-1]

for f in files:
    df = pd.read_csv('../results/'+f, sep='\t', index_col=0)
    df.index.name = 'gene_id'
    df = df.fillna(value=0)
    df = df.drop_duplicates()
    df.to_csv('../results/'+f, sep='\t')