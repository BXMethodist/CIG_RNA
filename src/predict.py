import pandas as pd, numpy as np
from scipy.stats import norm
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression

def label_label(table, positive_table):
    result = [0]*table.shape[0]
    table['label'] = result
    table.loc[table.index.isin(positive_table.index), 'label'] = 1
    return table

def center_normalization(table):
    return preprocessing.StandardScaler().fit(table)

def preprocessing_table(scaler, table):
    new_table_data = scaler.transform(table)
    new_table = pd.DataFrame(new_table_data, index=table.index, columns=table.columns)
    return new_table

def predict_LogisticRegression(X_train, Y_train, C=1., penalty='l1'):
    """
    :param table: a table which last columns is the label, (0 or 1), containing cell identity genes and non-cell identity genes
    :return: the sklearn logistic regression object
    """
    # print X_train
    predictor = LogisticRegression(penalty=penalty, C=C).fit(X_train, Y_train)
    return predictor

def predict_decision(predictor, real_table, label=False):
    if label:
        table = real_table.ix[:, :-1].copy()
    else:
        table = real_table.copy()
    result = predictor.decision_function(table)
    return pd.DataFrame(result, index=table.index)

def predict_proba(predictor, real_table, label=False):
    if label:
        table = real_table.ix[:, :-1].copy()
    else:
        table = real_table.copy()
    result = predictor.predict_proba(table)
    return pd.DataFrame(result, index=table.index)


if __name__ == '__main__':

    ## predict the lncRNA with CIG model and get the top list
    if False:
        CIG_df = pd.read_excel('../train/CIG_Bo_curated_311.xlsx', index_col=0)
        train_df = pd.read_csv('../train/CIG_grid_training_v311_size500-_v100.csv', index_col=0)
        columns = [col for col in train_df.columns if col.find('h3k')!=-1 and col.find('single')==-1 and col.find('genebody')==-1]
        columns1 = [col for col in columns if (col.find('genebody')==-1 or (col.find('coverage')!=-1)) and (col.find('h3k4me3')!=-1 or col.find('h3k27ac')!=-1)]
        columns2 = [col for col in columns if (col.find('genebody') != -1 or (col.find('coverage') != -1)) and (
        col.find('h3k4me1') != -1)]
        columns = columns1 + columns2
        print columns
        train_df = train_df[columns]

        train_df = label_label(train_df, CIG_df)

        scaler = center_normalization(train_df.iloc[:, :-1])

        train_df.iloc[:, :-1] = preprocessing_table(scaler, train_df.iloc[:, :-1])
        predictor = predict_LogisticRegression(train_df.iloc[:, :-1], train_df.iloc[:, -1], C=.2)

        print predictor.coef_
        print columns

        cell_types = ['HUVEC', 'CD34', 'GM12878', 'H1-hESC', 'HMEC', 'HSMM', 'MRG', 'MSC', 'neural', 'NHLF', ]
        for c in cell_types:
            real_table = pd.read_csv('../matrix/'+ c + '_parameters_real_table.csv', dtype={'gene_id': str})
            real_table = real_table.set_index('gene_id')
            real_table = real_table[columns]
            real_table = preprocessing_table(scaler, real_table)
            result = predict_decision(predictor, real_table)
            result2 = predict_proba(predictor, real_table)
            result = np.concatenate([result, result2], axis=1)
            result_df = pd.DataFrame(result, index=real_table.index)
            result_df.columns = ['distance', 'non-CIG_prob', 'CIG_prob']
            del result_df['non-CIG_prob']

            stats = importr('stats')
            result_df["p_value"] = 1 - norm.cdf(result_df['distance'], scale=1.25)
            result_df['FDR'] = stats.p_adjust(FloatVector(result_df["p_value"].tolist()),
                                       method='BH')

            result_df = result_df.sort_values(by=['distance'], ascending=False)
            result_df['rank'] = range(1, result_df.shape[0] + 1)
            result_df.to_csv('../results/'+ c + '_lncRNA_prediction.csv')

    ## Get the predicted candidates coordinate to view them in Genome browser
    if False:
        gtf = pd.read_csv('../ref_data/mitranscriptome_lncRNA.gtf', sep='\t', index_col=0)
        gtf['browser_index'] = gtf['hg19.knownGene.chrom']+':'+gtf['hg19.knownGene.txStart'].map(str)+'-'+gtf['hg19.knownGene.txEnd'].map(str)
        gtf['length'] = gtf['hg19.knownGene.txEnd']-gtf['hg19.knownGene.txStart']

        cell_types = ['HUVEC', 'CD34', 'GM12878', 'H1-hESC', 'HMEC', 'HSMM', 'MRG', 'MSC', 'neural', 'NHLF', ]
        for c in cell_types:
            df = pd.read_csv('../results/'+ c + '_lncRNA_prediction.csv', index_col=0)
            df['browser_index'] = gtf.ix[df.index, 'browser_index']
            df['length'] = gtf.ix[df.index, 'length']
            df.to_csv('../results/'+ c + '_lncRNA_prediction.csv')

    ## Get the predicted candidates with no overlapping of genes
    if True:
        lncRNA_no_overlap = pd.read_csv('../ref_data/mitranscriptome_lncRNA_no_overlap.bed', sep='\t', header=None)
        lncRNA_closest = pd.read_csv('../ref_data/mitranscriptome_lncRNA_closest.bed', sep='\t', header=None)
        lncRNA_closest.index = lncRNA_closest.iloc[:, 3]
        lncRNAs = lncRNA_no_overlap.iloc[:, 3].unique()
        cell_types = ['HUVEC', 'CD34', 'GM12878', 'H1-hESC', 'HMEC', 'HSMM', 'MRG', 'MSC', 'neural', 'NHLF', ]
        for c in cell_types:
            df = pd.read_csv('../results/' + c + '_lncRNA_prediction.csv', index_col=0)
            del df['not_overlap_with_genes']
            # del df['overlap_with_genes']
            df['not_overlap_with_genes'] = [0] * df.shape[0]
            df.ix[df.index.isin(lncRNAs), 'not_overlap_with_genes'] = 1
            lncRNA_closest.columns = ['lncRNA_chr', 'lncRNA_start', 'lncRNA_end', 'lncRNA_gene_id','lncRNA_dot', 'lncRNA_strand','closest_chr', 'closest_start', 'closest_end', 'closest_gene_id','closest_dot', 'closest_strand', 'closest_distance']
            df[['closest_chr', 'closest_start', 'closest_end', 'closest_gene_id','closest_dot', 'closest_strand', 'closest_distance']] = lncRNA_closest.ix[df.index, ['closest_chr', 'closest_start', 'closest_end', 'closest_gene_id','closest_dot', 'closest_strand', 'closest_distance']]
            del df['closest_dot']
            df = df[df['FDR']<0.01]
            df = df.sort_values(by=['not_overlap_with_genes', 'distance'], ascending=False)

            df.to_csv('../results/' + c + '_lncRNA_prediction.csv')
