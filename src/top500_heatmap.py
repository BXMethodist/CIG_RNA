import pandas as pd, numpy as np, os
from Wig import Wig
"""
This module is used for extracting the lncRNA signaling from the 10 different cell types wig. and make the heatmap
"""
def process_signal(signals):
    """
    change the wigChrom get signals to single types of array
    :param signals:
    :return:
    """
    index = signals[0, :]
    index = index.astype(int)
    start = np.min(index)
    end = np.max(index)

    new_signals = np.zeros(int(end-start+1))
    new_signals[index-start] = signals[1, :]
    return new_signals


def fetch_signals(wig, gene_infos, upstream=3000, downstream=10000, genebody=False):
    """

    :param wig: Wig object
    :param gene_infos: list of gene information (chr, start, end, strand)
    :param upstream: upstream of distance to the TSS
    :param downstream: downstream of distance to the TSS
    :return: a dataframe containing the value for that region
    """
    results = []
    index = []
    for gene_info in gene_infos:
        gene_id, chr, start, end, strand = gene_info

        if strand == '+' and not genebody:
            cur_start = start - upstream
            cur_end = start + downstream
            signal = wig.genome[chr].get_signals(cur_start, cur_end)
            signal = process_signal(signal)

        elif strand == '-' and not genebody:
            cur_start = end - downstream
            cur_end = end + upstream
            signal = wig.genome[chr].get_signals(cur_start, cur_end)
            signal = process_signal(signal)
            signal = np.flip(signal, axis=0)
        elif strand == '+' and genebody:
            cur_start = start - upstream
            cur_end = end
            signal = wig.genome[chr].get_signals(cur_start, cur_end)
            signal = process_signal(signal)
        elif strand == '-' and genebody:
            cur_start = start
            cur_end = end + upstream
            signal = wig.genome[chr].get_signals(cur_start, cur_end)
            signal = process_signal(signal)
            signal = np.flip(signal, axis=0)
        else:
            continue
        results.append(list(signal))
        index.append(gene_id)
    df = pd.DataFrame(results, index=index)
    return df

if __name__ == "__main__":

    ## To get the top 500 list for lncRNA with h3k4me3 or h3k27me3
    if False:
        files = ['../matrix/'+x for x in os.listdir('../matrix/') if x.endswith('real_table.xls')]

        gtf = pd.read_csv('../ref_data/mitranscriptome_lncRNA.gtf', sep='\t')
        gtf.index = gtf['hg19.knownGene.name']

        for f in files:
            cell_type = f.split('_')[0].replace('../matrix/', '')
            df = pd.read_csv(f, sep='\t', index_col=0)
            df.index.name = 'id'
            k4m3 = df.nlargest(500, 'h3k4me3_qn_total_width').name.values
            k27m3 = df.nlargest(500, 'h3k27me3_qn_total_width_genebody').name.values

            k4m3_df = gtf.ix[k4m3, ['hg19.knownGene.name', 'hg19.knownGene.chrom', 'hg19.knownGene.txStart',
                                    'hg19.knownGene.txEnd', 'hg19.knownGene.strand']].values

            k27m3_df = gtf.ix[k27m3, ['hg19.knownGene.name', 'hg19.knownGene.chrom', 'hg19.knownGene.txStart',
                                    'hg19.knownGene.txEnd', 'hg19.knownGene.strand']].values

            cur_k4m3_wig = Wig('/scratch/tmhbxx3/CIG/h3k4me3_qn/'+cell_type+'_h3k4me3.qnor.wig')
            cur_k27m3_wig = Wig('/scratch/tmhbxx3/CIG/h3k27me3_qn/' + cell_type + '_h3k27me3.qnor.wig')

            # get the signal for broad h3k4me3 lncRNA
            cur_k4m3_result = fetch_signals(cur_k4m3_wig, k4m3_df)
            cur_k27m3_result = fetch_signals(cur_k27m3_wig, k4m3_df)

            cur_k4m3_result.to_csv('../results/'+cell_type+'_h3k4me3_lncRNA_k4m3sig.xls', sep='\t')
            cur_k27m3_result.to_csv('../results/'+cell_type + '_h3k4me3_lncRNA_k27m3sig.xls', sep='\t')

            # get the signal for broad h3k27me3 lncRNA
            cur_k4m3_result = fetch_signals(cur_k4m3_wig, k27m3_df)
            cur_k27m3_result = fetch_signals(cur_k27m3_wig, k27m3_df)

            cur_k4m3_result.to_csv('../results/'+cell_type + '_h3k27me3_lncRNA_k4m3sig.xls', sep='\t')
            cur_k27m3_result.to_csv('../results/' + cell_type + '_h3k27me3_lncRNA_k27m3sig.xls', sep='\t')

    ## To get the top 500 list for genes with h3k4me3 or h3k27me3
    if False:
        files = ['../matrix_gene/' + x for x in os.listdir('../matrix_gene/') if x.endswith('real_table.xls')]

        gtf = pd.read_csv('../ref_data/hg19.GREATgene2UCSCknownGenes.table.xls', sep='\t')
        gtf.index = gtf['hg19.kgXref.geneSymbol']

        for f in files:
            cell_type = f.replace('../matrix_gene/', '').split('_')[0]
            df = pd.read_csv(f, sep='\t', index_col=0)
            df.index.name = 'id'

            df = df[df.name.isin(gtf.index)]

            k4m3 = df.nlargest(500, 'h3k4me3_qn_total_width').name.values
            k27m3 = df.nlargest(500, 'h3k27me3_qn_total_width_genebody').name.values

            k4m3_df = gtf.ix[k4m3, ['hg19.kgXref.geneSymbol', 'hg19.knownGene.chrom', 'hg19.knownGene.txStart',
                                    'hg19.knownGene.txEnd', 'hg19.knownGene.strand']].values

            k27m3_df = gtf.ix[k27m3, ['hg19.knownGene.name', 'hg19.knownGene.chrom', 'hg19.knownGene.txStart',
                                      'hg19.knownGene.txEnd', 'hg19.knownGene.strand']].values

            cur_k4m3_wig = Wig('/scratch/tmhbxx3/CIG/h3k4me3_qn/' + cell_type + '_h3k4me3.qnor.wig')
            cur_k27m3_wig = Wig('/scratch/tmhbxx3/CIG/h3k27me3_qn/' + cell_type + '_h3k27me3.qnor.wig')

            # get the signal for broad h3k4me3 lncRNA
            cur_k4m3_result = fetch_signals(cur_k4m3_wig, k4m3_df)
            cur_k27m3_result = fetch_signals(cur_k27m3_wig, k4m3_df)

            cur_k4m3_result.to_csv('../results/' + cell_type + '_h3k4me3_gene_k4m3sig.xls', sep='\t')
            cur_k27m3_result.to_csv('../results/' + cell_type + '_h3k4me3_gene_k27m3sig.xls', sep='\t')

            # get the signal for broad h3k27me3 lncRNA
            cur_k4m3_result = fetch_signals(cur_k4m3_wig, k27m3_df)
            cur_k27m3_result = fetch_signals(cur_k27m3_wig, k27m3_df)

            cur_k4m3_result.to_csv('../results/' + cell_type + '_h3k27me3_gene_k4m3sig.xls', sep='\t')
            cur_k27m3_result.to_csv('../results/' + cell_type + '_h3k27me3_gene_k27m3sig.xls', sep='\t')



