"""
check whether CSE and OSE are enriched in cancer pathway
"""

import pandas as pd, numpy as np

def write_txt(list, name):
    f = open(name, 'w')
    for l in list:
        f.write(l+'\n')
    f.close()


if __name__ == "__main__":
    if True:
        celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte',
                     'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
                     'neural-cell', 'osteoblast', 'smooth-muscle-cell']
        # 'H1-hESC',
        # gtf_df = pd.read_csv('/home/tmhbxx3/archive/ref_data/hg19/hg19.GREATgene2UCSCknownGenes.table.xls',
        #                      sep='\t', index_col=0)
        gtf_df = pd.read_csv('../ref_data/hg19.GREATgene2UCSCknownGenes.table.xls',
                             sep='\t', index_col=0)

        ogs = pd.read_csv('../ref_data/OG_top500.xls', sep='\t', index_col=0).index
        tsgs = pd.read_csv('../ref_data/TSG_top500.xls', sep='\t', index_col=0).index

        cse_df = pd.DataFrame(index=gtf_df.index)
        ose_df = pd.DataFrame(index=gtf_df.index)
        control_df = pd.DataFrame(index=gtf_df.index)

        for celltype in celltypes:
            print celltype
            try:
                cse_genes = [x.strip() for x in set(pd.read_csv('../CSE_genes/'+celltype + '_CSE_genes.xls', sep='\t')['gene_id'].unique())]
                ose_genes = [x.strip() for x in set(pd.read_csv('../CSE_genes/'+celltype + '_OSE_genes.xls', sep='\t')['gene_id'].unique())]
            except:
                cse_genes = [x.strip() for x in set(pd.read_csv('../CSE_genes/'+celltype + '_CSE_genes.xls', sep='\t', header=None).iloc[:, 3].unique())]
                ose_genes = [x.strip() for x in set(pd.read_csv('../CSE_genes/'+celltype + '_OSE_genes.xls', sep='\t', header=None).iloc[:, 3].unique())]
            cur_control_genes = set(gtf_df.sample((len(cse_genes) + len(ose_genes)) / 2).index)

            cse_df.ix[list(cse_genes), celltype] = 1
            ose_df.ix[list(ose_genes), celltype] = 1
            control_df.ix[list(cur_control_genes), celltype] =1

            print celltype
            print len(set(cse_genes).intersection(tsgs))
            print len(set(ose_genes).intersection(tsgs))

            print len(set(cse_genes).intersection(ogs))
            print len(set(ose_genes).intersection(ogs))

        conserved_cse = cse_df[cse_df.sum(axis=1)>10].index
        unique_cse = cse_df[cse_df.sum(axis=1)<5].index

        conserved_ose = ose_df[ose_df.sum(axis=1) > 10].index
        unique_ose = ose_df[ose_df.sum(axis=1) < 5].index

        conserved_control = control_df[control_df.sum(axis=1) > 10].index
        unique_control = control_df[control_df.sum(axis=1) < 5].index

        write_txt(conserved_cse, '../CSE_genes/conserved_cse.txt')
        write_txt(conserved_ose, '../CSE_genes/conserved_ose.txt')
        write_txt(conserved_control, '../CSE_genes/conserved_control.txt')
        write_txt(unique_cse, '../CSE_genes/unique_cse.txt')
        write_txt(unique_ose, '../CSE_genes/unique_ose.txt')
        write_txt(unique_control, '../CSE_genes/unique_control.txt')

        print len(conserved_cse.intersection(tsgs))
        print len(conserved_ose.intersection(tsgs))

        print len(conserved_cse.intersection(ogs))
        print len(conserved_ose.intersection(ogs))

        print len(unique_cse.intersection(tsgs))
        print len(unique_ose.intersection(tsgs))

        print len(unique_cse.intersection(ogs))
        print len(unique_ose.intersection(ogs))




