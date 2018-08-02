import pandas as pd, os

if __name__ == '__main__':
    gene_df = pd.read_csv('../ref_data/hg19.ucscgenes.knowngene_proteinOnly.xls', sep='\t', index_col=0)
    lncRNA_df = pd.read_csv('../ref_data/mitranscriptome_lncRNA.gtf', sep='\t', index_col=0)
    gene_df['name'] = gene_df['hg19.kgXref.geneSymbol']
    gene_df['score'] = ['.']* gene_df.shape[0]
    lncRNA_df['name'] = lncRNA_df.index
    lncRNA_df['score'] = ['.'] * lncRNA_df.shape[0]

    gene_df[[u'hg19.knownGene.chrom', u'hg19.knownGene.txStart', u'hg19.knownGene.txEnd', 'name', 'score', 'hg19.knownGene.strand']].to_csv('../ref_data/hg19.ucscgenes.knowngene_proteinOnly.bed',sep='\t', header=None, index = None)
    lncRNA_df[[u'hg19.knownGene.chrom', u'hg19.knownGene.txStart', u'hg19.knownGene.txEnd', 'name', 'score', 'hg19.knownGene.strand']].to_csv('../ref_data/mitranscriptome_lncRNA.bed',sep='\t', header=None, index = None)

    os.system('/Users/boxia/tools/bedtools2/bin/bedtools intersect -v -a mitranscriptome_lncRNA.bed -b hg19.ucscgenes.knowngene_proteinOnly.bed > mitranscriptome_lncRNA_no_overlap.bed')

    os.system('sort -k1,1 -k2,2n hg19.ucscgenes.knowngene_proteinOnly.bed > hg19.ucscgenes.knowngene_proteinOnly.sorted.bed')
    os.system('/Users/boxia/tools/bedtools2/bin/bedtools closest -t first -d -a mitranscriptome_lncRNA.bed -b hg19.ucscgenes.knowngene_proteinOnly.sorted.bed > mitranscriptome_lncRNA_closest.bed')