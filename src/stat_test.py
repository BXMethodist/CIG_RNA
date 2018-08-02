from scipy import stats
import pandas as pd, numpy as np

def logP_wilcoxon(groupA, groupB, bags=1):
    """
    return the negative log P value for two groups
    :param groupA:
    :param groupB:
    :return:
    """
    p_values = []

    for i in range(bags):
        # cur_groupA = np.random.choice(groupA, int(len(groupA)*0.75), replace=True)
        # cur_groupB = np.random.choice(groupB, int(len(groupB)*0.75), replace=True)
        cur_groupA = groupA
        cur_groupB = groupB
        try:
            rank_diff, p = stats.mannwhitneyu(cur_groupA, cur_groupB, alternative='less')
            # print p
            # There is a problem in multiple processing. for some reason, the function will return a very small number
            # However, it could pass the test without problem if just run the corresponding dataframe
            # if np.log10(p) < -80:
            #     p_values.append(0)

            p_values.append(np.log10(p))
        except:
            p_values.append(0)

    return np.mean(p_values)

def logP_fisher(gene_df, all_stat_df, criteria, top_enrich=500, ascending=False):
    total_genes = all_stat_df.shape[0]
    sort_result = all_stat_df.sort_values(by=[criteria], ascending=ascending)
    top_genes = sort_result['gene'].tolist()[:top_enrich]
    overlap = len(set(top_genes).intersection(set(gene_df['gene'].tolist())))
    not_overlap = top_enrich - overlap
    # print overlap, top_enrich, not_overlap, total_genes
    p = stats.fisher_exact([[overlap, not_overlap], [not_overlap, total_genes - 2 * top_enrich + overlap]],
                           alternative='greater')[1]
    return np.log10(p)


if __name__ == "__main__":
    total_genes = 28000
    overlap = 43
    not_overlap =457
    top_enrich = 500

    p = stats.fisher_exact([[overlap, not_overlap], [not_overlap, total_genes - 2 * top_enrich + overlap]],
                           alternative='greater')[1]

    print -np.log10(p)
    #
    # p = stats.fisher_exact([[1, 4], [3, 12 - 2 * 4 + 1]],
    #                        alternative='greater')[1]
    # print p

    print stats.fisher_exact([[14, 72], [55, 230]], alternative='less')[1]