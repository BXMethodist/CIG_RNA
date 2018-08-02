"""
CTCF is the one that can provide additional information
"""

import pandas as pd, numpy as np, os

if __name__ == "__main__":
    """
    whether CTCF can provide an additional layer of information by using correlation
    """
    if True:
        final_df = pd.DataFrame()
        tables = ['../real_tables/'+x for x in os.listdir('../real_tables/') if x.endswith('.csv')]
        for table in tables:
            celltype = table[table.rfind('/')+1:].split('_')[0]
            df = pd.read_csv(table, index_col=0)
            columns = [c for c in df.columns if c.find('height')!=-1 or c.find('total_width')!=-1]
            df = df[columns]
            df = df[df['H3K27ac_total_width']>=500]
            corr = df.corr(method='spearman')
            corr = corr.ix[:, 'H3K27ac_total_width']
            final_df[celltype] = corr
        final_df['avg'] = final_df.mean(axis=1)
        final_df.to_csv('../CTCF/correlation_H3K27ac_total_width.xls', sep='\t')

