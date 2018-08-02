import pandas as pd, os, numpy as np

if __name__ == "__main__":
    if True:
        celltypes = ['HSMM', 'HMEC', 'HUVEC',  'NHLF']

        """
        'cardiac-muscle','H1-hESC','neutrophil','CD4', 'B-cell','monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
        'osteoblast', 'smooth-muscle-cell',
        """

        fpkm = pd.read_csv('../RNAexp/genes.fpkm.txt', sep='\t', index_col=0)

        for celltype in celltypes:
            CSE_genes = pd.read_csv('../CSE_genes/'+celltype+'_CSE_genes.xls', sep='\t')['gene_id'].unique()
            OSE_genes = pd.read_csv('../CSE_genes/'+celltype+'_OSE_genes.xls', sep='\t')['gene_id'].unique()

            cse_stem_df = fpkm[fpkm.index.isin(CSE_genes)][['H1-hESC_FPKM', celltype+'_FPKM']]
            ose_stem_df = fpkm[fpkm.index.isin(OSE_genes)][['H1-hESC_FPKM', celltype+'_FPKM']]

            cse_stem_df.to_csv(celltype+'_stem_CSE.xls', sep='\t')
            ose_stem_df.to_csv(celltype + '_stem_OSE.xls', sep='\t')

            length = max(len(CSE_genes), len(OSE_genes))
            final_df = pd.DataFrame(index=range(max(len(CSE_genes), len(OSE_genes))))
            print final_df.shape[0]
            print len(list(fpkm[fpkm.index.isin(CSE_genes)][celltype+'_FPKM'].values))

            final_df['CSE'] = list(fpkm[fpkm.index.isin(CSE_genes)][celltype+'_FPKM'].values) + [np.nan]* (length- len(list(fpkm[fpkm.index.isin(CSE_genes)][celltype+'_FPKM'].values)))
            final_df['OSE'] = list(fpkm[fpkm.index.isin(OSE_genes)][celltype+'_FPKM'].values) + [np.nan]* (length- len(list(fpkm[fpkm.index.isin(OSE_genes)][celltype+'_FPKM'].values)))
            final_df.to_csv(celltype+'_CSE_vs_OSE.xls', sep='\t')


