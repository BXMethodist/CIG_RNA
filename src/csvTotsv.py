import pandas as pd, numpy as np, os

def csvTotsv(csv, outpath):
    cell_type = csv.split('_')[0]
    df = pd.read_csv(outpath+csv, index_col=0)
    # index = df.index + '_' + cell_type
    # print index
    df['name'] = df.index
    df.index = df.index +'_'+cell_type
    df.to_csv(outpath+csv[:-4]+'.xls', sep='\t')
    return df


if __name__ == "__main__":
    ## for lncRNA
    if False:
        dfs = []
        files = [x for x in os.listdir('../matrix/') if x.endswith('.csv')]
        for f in files:
            print f
            cur_df = csvTotsv(f, '../matrix/')
            dfs.append(cur_df)
        final_df = dfs[0]
        for i in range(1, len(dfs)):
            final_df = final_df.append(dfs[i])
        final_df.to_csv('ten_cell_types_lncRNA.xls', sep='\t')

    if False:
        dfs = []
        files = [x for x in os.listdir('../matrix_gene/') if x.endswith('.csv')]
        for f in files:
            cur_df = csvTotsv(f, '../matrix_gene/')
            dfs.append(cur_df)
        final_df = dfs[0]
        for i in range(1, len(dfs)):
            final_df = final_df.append(dfs[i])
        final_df.to_csv('ten_cell_types_lncRNA_gene.xls', sep='\t')
