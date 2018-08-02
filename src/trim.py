import pandas as pd

## Extract the lncRNA annotation from gtf
if __name__ == "__main__":
    df = pd.read_csv('../ref_data/mitranscriptome.v2.gtf', sep='\t', header=None)

    df = df[~((df.iloc[:, 8].str.contains('pseudogene'))|(df.iloc[:, 8].str.contains('exon_number'))|(df.iloc[:, 8].str.contains('tucp'))|(df.iloc[:, 8].str.contains('mixed'))|(df.iloc[:, 8].str.contains('protein_coding')))]

    df.to_csv('../ref_data/mitranscriptome_nopseudo_exon.gtf', index=None, header=None, sep='\t')

    names = []

    for i in range(df.shape[0]):
        name = df.iloc[i, 8].split(';')
        for n in name:
            if n.find('transcript_id')!=-1:
                n = n.split()[-1].strip().replace('"', '')
                names.append(n)
                break

    df.index = names
    df.index.name = 'hg19.knownGene.name'
    df.columns = ['hg19.knownGene.chrom', 'name', 'category', 'hg19.knownGene.txStart', 'hg19.knownGene.txEnd',  'value',  'hg19.knownGene.strand', 'dot', 'information']
    df['hg19.kgXref.geneSymbol'] = list(df.index)


    df.to_csv('../ref_data/mitranscriptome_lncRNA.gtf', sep='\t')