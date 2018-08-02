"""
check whether CSE have more TFs compared to OSE
"""
import pandas as pd, numpy as np, os
from scipy import stats

if __name__ == "__main__":
    """
    compare the number of TFs
    """
    if False:
        TFs = set(pd.read_csv('transcript_factors.txt', sep='\t', index_col=0).index)
        celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte',
                     'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
                     'neural-cell', 'osteoblast', 'smooth-muscle-cell']
        #'H1-hESC',
        gtf_df = pd.read_csv('/home/tmhbxx3/archive/ref_data/hg19/hg19.GREATgene2UCSCknownGenes.table.xls',
                             sep='\t', index_col=0)

        results = []
        for celltype in celltypes:
            print celltype
            try:
                cse_genes = set(pd.read_csv(celltype + '_CSE_genes.xls', sep='\t')['gene_id'].unique())
                ose_genes = set(pd.read_csv(celltype+'_OSE_genes.xls', sep='\t')['gene_id'].unique())
            except:
                cse_genes = set(pd.read_csv(celltype + '_CSE_genes.xls', sep='\t', header=None).iloc[:,3].unique())
                ose_genes = set(pd.read_csv(celltype + '_OSE_genes.xls', sep='\t', header=None).iloc[:,3].unique())
            cur_control_genes = set(gtf_df.sample((len(cse_genes)+len(ose_genes))/2).index)

            cse = len(TFs.intersection(cse_genes))
            ose = len(TFs.intersection(ose_genes))
            control = len(TFs.intersection(cur_control_genes))
            print celltype, cse, ose, control
            results.append((celltype, cse, ose, control))
        result_df = pd.DataFrame(results)
        result_df.columns = ['celltype', 'CSE', 'OSE', 'control']
        result_df.to_csv('TFs_CSE_vs_OSE.xls', sep='\t', index=None)
    """
        compare the number of TFs, check p value
    """
    if False:
        df = pd.read_csv('../CSE_genes/TFs_CSE_vs_OSE.xls', sep='\t', index_col=0)
        print stats.mannwhitneyu(df['control'], df['OSE'], alternative='greater')[1]
        print stats.mannwhitneyu(df['control'], df['CSE'], alternative='less')[1]

    """
    check whether more network edges:
    """
    if False:
        GRN_df = pd.read_excel('Human_Big_GRN_032014.xlsx')
        GRN_df = GRN_df[GRN_df['type'] == 'all']

        celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte',
                     'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
                     'neural-cell', 'osteoblast', 'smooth-muscle-cell']
        gtf_df = pd.read_csv('/home/tmhbxx3/archive/ref_data/hg19/hg19.GREATgene2UCSCknownGenes.table.xls',
                             sep='\t', index_col=0)

        results = []
        for celltype in celltypes:
            print celltype
            try:
                cse_genes = set(pd.read_csv(celltype + '_CSE_genes.xls', sep='\t')['gene_id'].unique())
                ose_genes = set(pd.read_csv(celltype + '_OSE_genes.xls', sep='\t')['gene_id'].unique())
            except:
                cse_genes = set(pd.read_csv(celltype + '_CSE_genes.xls', sep='\t', header=None).iloc[:, 3].unique())
                ose_genes = set(pd.read_csv(celltype + '_OSE_genes.xls', sep='\t', header=None).iloc[:, 3].unique())
            cur_control_genes = set(gtf_df.sample((len(cse_genes) + len(ose_genes)) / 2).index)


            cse_edges = GRN_df[(GRN_df['TG'].isin(cse_genes)) & (GRN_df['TF'].isin(cse_genes))].shape[0]
            ose_edges = GRN_df[(GRN_df['TG'].isin(ose_genes)) & (GRN_df['TF'].isin(ose_genes))].shape[0]
            control_edges = GRN_df[(GRN_df['TG'].isin(cur_control_genes)) & (GRN_df['TF'].isin(cur_control_genes))].shape[0]

            results.append((celltype, cse_edges, ose_edges, control_edges))
        result_df = pd.DataFrame(results)
        result_df.columns = ['celltype', 'CSE', 'OSE', 'control']
        result_df.to_csv('network_edgs_CSE_OSE.xls', sep='\t', index=None)

    if False:
        df = pd.read_csv('../CSE_genes/network_edgs_CSE_OSE.xls', sep='\t', index_col=0)
        print stats.mannwhitneyu(df['control'], df['OSE'], alternative='less')[1]
        print stats.mannwhitneyu(df['control'], df['CSE'], alternative='less')[1]








