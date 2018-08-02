"""
define integrated activation domain, defined by h3k4me1, h3k4me3, h3k27ac broad domain, CTCF
"""
import pandas as pd, os, numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.preprocessing import MinMaxScaler

def tau(final_df):
    cdf = final_df + abs(final_df.min().min())
    return (1 - cdf.divide(cdf.max(axis=1), axis=0)).sum(axis=1) / (cdf.shape[1] - 1)

def merge(df):
    results = []
    cur_chr, cur_start, cur_end, cur_center, cur_width_above_cutoff, cur_total_signal, cur_height, cur_height_logP = \
        None, None, None, None, None, None, None, None
    for i in range(df.shape[0]):
        chr, start, end, center, width_above_cutoff, total_signal, height, height_logP = df.ix[i, :]
        if cur_chr is None:
            cur_chr, cur_start, cur_end, cur_center, cur_width_above_cutoff, cur_total_signal, cur_height, cur_height_logP = chr, start, end, center, width_above_cutoff, total_signal, height, height_logP
        elif cur_chr != chr:
            results.append([cur_chr, cur_start, cur_end, cur_center, cur_width_above_cutoff, cur_total_signal, cur_height, cur_height_logP])
            cur_chr, cur_start, cur_end, cur_center, cur_width_above_cutoff, cur_total_signal, cur_height, cur_height_logP = chr, start, end, center, width_above_cutoff, total_signal, height, height_logP
        elif start - cur_end > 3000:
            results.append(
                [cur_chr, cur_start, cur_end, cur_center, cur_width_above_cutoff, cur_total_signal, cur_height,
                 cur_height_logP])
            cur_chr, cur_start, cur_end, cur_center, cur_width_above_cutoff, cur_total_signal, cur_height, cur_height_logP = chr, start, end, center, width_above_cutoff, total_signal, height, height_logP
        else:
            cur_end = end
            cur_center = (cur_start + cur_end)/2
            cur_width_above_cutoff += width_above_cutoff
            cur_total_signal += total_signal
            cur_height = cur_height if height < cur_height else height
            cur_height_logP = cur_height_logP if height_logP < cur_height_logP else height_logP
    results.append([cur_chr, cur_start, cur_end, cur_center, cur_width_above_cutoff, cur_total_signal, cur_height,
                 cur_height_logP])
    final_df = pd.DataFrame(results)
    final_df.columns = df.columns
    return final_df

def SE_cutoff(df):
    cur_df = df.copy()
    cur_df = cur_df.sort_values(by='SE_width', ascending=True)
    cur_df['rank'] = range(cur_df.shape[0])
    scaler = MinMaxScaler()
    scaler.fit(cur_df['SE_width'].values.reshape(-1, 1))
    cur_df['transformed_width'] = scaler.transform(cur_df['SE_width'].values.reshape(-1, 1))

    scaler = MinMaxScaler()
    scaler.fit(cur_df['rank'].values.reshape(-1, 1))
    cur_df['transformed_rank'] = scaler.transform(cur_df['rank'].values.reshape(-1, 1))

    for i in range(cur_df.shape[0]):
        a, b = 1, cur_df.ix[i, 'transformed_width'] - cur_df.ix[i, 'transformed_rank']
        cur_df['new_width'] = cur_df['transformed_rank'] + b
        if cur_df[cur_df['new_width'] < cur_df['transformed_width']].shape[0] == 0:
            return i
    return

def random_df(df, chr_sizes):
    results = []
    columns = df.columns
    for i in range(df.shape[0]):
        try:
            chr, start, end = df.ix[i, 'chr'], df.ix[i, 'start'], df.ix[i, 'end']
            center, width, total_signal, height, height_logP = df.ix[i, 'center'], \
                                                               df.ix[i, 'width_above_cutoff'], \
                                                               df.ix[i, 'total_signal'], \
                                                               df.ix[i, 'height'], \
                                                               df.ix[i, 'height_logP'], \

            size = chr_sizes.ix[chr, 'size']
            size = size - (end-start)
            new_start = np.random.choice(size)
            results.append([chr, new_start, new_start+(end-start), center, width, total_signal, height, height_logP])
        except:
            continue
    result_df = pd.DataFrame(results)
    result_df.columns = columns
    return result_df


if __name__ == "__main__":
    """
    How much CTCF peak overlap with super enhancer
    """
    if False:
        os.system('module load bedtools')
        # celltypes = os.listdir('./peaks/')
        celltypes = ['H1-hESC']
        for celltype in celltypes:
            # 1 remove all CTCF peak below 100
            # 2 merge H3K27ac peaks within 3kb
            # 3 overlap CTCF and H3K27ac with different distance gap
            CTCF_df = pd.read_csv('./peaks/'+celltype+'/CTCF/'+celltype+'_CTCF_10.xls', sep='\t')
            CTCF_df = CTCF_df[CTCF_df['height'] >= 100]
            CTCF_df.to_csv('./peaks/'+celltype+'/CTCF/'+celltype+'_CTCF_10_heigh100.xls', sep='\t', index=None, header=None)

            os.system(
                'sort -k1,1 -k2,2n ./peaks/'+celltype+'/CTCF/'+celltype+'_CTCF_10_heigh100.xls > ./peaks/'+celltype+'/CTCF/'+celltype+'_CTCF_10_sorted_heigh100.xls')

            os.system('sort -k1,1 -k2,2n ./peaks/'+celltype+'/H3K27ac/'+celltype+'_H3K27ac_10.xls > ./peaks/'+celltype+'/H3K27ac/'+celltype+'_H3K27ac_10_sorted.xls')

            SE_df = pd.read_csv('./peaks/'+celltype+'/H3K27ac/'+celltype+'_H3K27ac_10_sorted.xls', sep='\t')

            SE_df = merge(SE_df)

            SE_df.to_csv('./peaks/'+celltype+'/H3K27ac/'+celltype+'_H3K27ac_10_merge.xls', sep='\t', index=None, header=None)

            os.system('bedtools closest -a ./peaks/'+celltype+'/H3K27ac/'+celltype+'_H3K27ac_10_merge.xls -b ./peaks/'+celltype+'/CTCF/'+celltype+'_CTCF_10_sorted_heigh100.xls > '+celltype+'_CTCF_H3K27ac_refmap.xls' )

            CTCF_SE_df = pd.read_csv(celltype+'_CTCF_H3K27ac_refmap.xls', sep='\t', header=None)
            results = []
            for k in range(CTCF_SE_df.shape[0]):
                SE_chr, SE_start, SE_end, CTCF_start, CTCF_end, CTCF_height = CTCF_SE_df.iloc[k, [0, 1, 2, 9, 10, 14]]
                SE_width = SE_end - SE_start
                if SE_start <= CTCF_start <= SE_end or SE_start <= CTCF_end <= SE_end:
                    distance = 0
                else:
                    distance = abs(SE_end - CTCF_start) if abs(SE_end - CTCF_start) < abs(SE_start - CTCF_end) \
                        else abs(SE_start - CTCF_end)
                results.append([SE_chr, SE_start, SE_end, CTCF_start, CTCF_end, SE_width, CTCF_height, distance])
            final_df = pd.DataFrame(results)
            final_df.columns = ['SE_chr', 'SE_start', 'SE_end', 'CTCF_start', 'CTCF_end',
                                'SE_width', 'CTCF_height', 'distance']
            final_df.to_csv(celltype+'_CTCF_H3K27ac_refmap.xls', sep='\t', index=None)

    """
        How much CTCF peak overlap with random control super enhancer
    """
    if False:
        chr_sizes = pd.read_csv('/home/tmhbxx3/archive/ref_data/hg19/hg19_chr_sizes.txt', sep='\t', header=None)
        chr_sizes.columns =['chr', 'size']
        chr_sizes = chr_sizes.set_index(['chr'])
        os.system('module load bedtools')
        celltypes = os.listdir('./peaks/')
        # celltypes = ['H1-hESC']
        for celltype in celltypes:
            # 1 remove all CTCF peak below 100
            # 2 merge H3K27ac peaks within 3kb
            # 3 overlap CTCF and H3K27ac with different distance gap
            CTCF_df = pd.read_csv('./peaks/' + celltype + '/CTCF/' + celltype + '_CTCF_10.xls', sep='\t')
            CTCF_df = CTCF_df[CTCF_df['height'] >= 100]
            CTCF_df.to_csv('./peaks/' + celltype + '/CTCF/' + celltype + '_CTCF_10_heigh100.xls', sep='\t', index=None,
                           header=None)

            os.system(
                'sort -k1,1 -k2,2n ./peaks/' + celltype + '/CTCF/' + celltype + '_CTCF_10_heigh100.xls > ./peaks/' + celltype + '/CTCF/' + celltype + '_CTCF_10_sorted_heigh100.xls')

            se_df = pd.read_csv('./peaks/' + celltype + '/H3K27ac/' + celltype + '_H3K27ac_10.xls', sep='\t')
            se_df = random_df(se_df, chr_sizes)
            se_df.to_csv('./peaks/' + celltype + '/H3K27ac/' + celltype + '_H3K27ac_10_random.xls', sep='\t', index=None)

            os.system(
                'sort -k1,1 -k2,2n ./peaks/' + celltype + '/H3K27ac/' + celltype + '_H3K27ac_10_random.xls > ./peaks/' + celltype + '/H3K27ac/' + celltype + '_H3K27ac_10_random_sorted.xls')

            SE_df = pd.read_csv('./peaks/' + celltype + '/H3K27ac/' + celltype + '_H3K27ac_10_random_sorted.xls', sep='\t')

            SE_df = merge(SE_df)

            SE_df.to_csv('./peaks/' + celltype + '/H3K27ac/' + celltype + '_H3K27ac_10_random_merge.xls', sep='\t', index=None,
                         header=None)

            os.system(
                'bedtools closest -a ./peaks/' + celltype + '/H3K27ac/' + celltype + '_H3K27ac_10_random_merge.xls -b ./peaks/' + celltype + '/CTCF/' + celltype + '_CTCF_10_sorted_heigh100.xls > ' + celltype + '_CTCF_H3K27ac_refmap_random.xls')
            #bedtools closest -a ./peaks/H1-hESC/H3K27ac/H1-hESC_H3K27ac_10_random_merge.xls -b ./peaks/H1-hESC/CTCF/H1-hESC_CTCF_10_sorted_heigh100.xls > H1-hESC_CTCF_H3K27ac_refmap_random.xls'

            CTCF_SE_df = pd.read_csv(celltype + '_CTCF_H3K27ac_refmap_random.xls', sep='\t', header=None)
            results = []
            for k in range(CTCF_SE_df.shape[0]):
                SE_chr, SE_start, SE_end, CTCF_start, CTCF_end, CTCF_height = CTCF_SE_df.iloc[k, [0, 1, 2, 9, 10, 14]]
                SE_width = SE_end - SE_start
                if SE_start <= CTCF_start <= SE_end or SE_start <= CTCF_end <= SE_end:
                    distance = 0
                else:
                    distance = abs(SE_end - CTCF_start) if abs(SE_end - CTCF_start) < abs(SE_start - CTCF_end) \
                        else abs(SE_start - CTCF_end)
                results.append([SE_chr, SE_start, SE_end, CTCF_start, CTCF_end, SE_width, CTCF_height, distance])
            final_df = pd.DataFrame(results)
            final_df.columns = ['SE_chr', 'SE_start', 'SE_end', 'CTCF_start', 'CTCF_end',
                                'SE_width', 'CTCF_height', 'distance']
            final_df.to_csv(celltype + '_CTCF_H3K27ac_refmap_random.xls', sep='\t', index=None)
    """
    Get the random percentage
    """
    if True:
        files = ['../CTCF_SE/'+x for x in os.listdir('../CTCF_SE/') if x.endswith('_refmap_random.xls')]
        results = []
        for f in files:
            cur_df = pd.read_csv(f, sep='\t')
            te_p = cur_df[cur_df['distance']==0].shape[0]*1./cur_df.shape[0]
            cur_se_df = cur_df.nlargest(1000, 'SE_width')
            se_p = cur_se_df[cur_se_df['distance']==0].shape[0]*1./cur_se_df.shape[0]
            results.append((f.split('_')[0], te_p, se_p))
        result_df = pd.DataFrame(results)
        result_df.columns = ['celltype', 'TE/SE', 'SE']
        result_df.to_csv('../CTCF_SE/CTCF_H3K27ac_percentage_random.xls', sep='\t')

    """
    Only SE is associated with CTCF not TE
    """
    if False:
        celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte', 'cardiac-muscle',
                     'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
                     'neural-cell', 'neutrophil', 'osteoblast', 'smooth-muscle-cell', 'H1-hESC']
        results = []
        for celltype in celltypes:
            CTCF_SE_df = pd.read_csv('../CTCF_SE/'+celltype + '_CTCF_H3K27ac_refmap.xls', sep='\t')
            withCTCF = CTCF_SE_df[CTCF_SE_df['distance'] == 0].copy()
            woCTCF = CTCF_SE_df[CTCF_SE_df['distance'] != 0].copy()
            CTCF_SE_df['distance'] = np.log10(CTCF_SE_df['distance']+1)
            CTCF_SE_df['SE_width'] = np.log10(CTCF_SE_df['SE_width'])
            plot = CTCF_SE_df.plot.scatter(x='SE_width', y='distance', edgecolors=None, c='red', s=1)
            fig = plot.get_figure()
            fig.savefig('../CTCF_SE/'+celltype + '_CTCF_H3K27ac_width_vs_distance.pdf')
            plt.close("all")

            width_df = pd.DataFrame(index=range(max(withCTCF.shape[0], woCTCF.shape[0])), columns=['CTCF', 'none CTCF'])
            width_df['CTCF'] = list(withCTCF['SE_width'].values) + [np.nan]* (width_df.shape[0] - withCTCF.shape[0])
            width_df['none CTCF'] = list(woCTCF['SE_width'].values) + [np.nan] * (width_df.shape[0] - woCTCF.shape[0])

            width_df.to_csv('../CTCF_SE/'+celltype + '_CTCF_H3K27ac_width_cmp.xls', sep='\t', index=None)

            p_value = stats.mannwhitneyu(width_df['CTCF'], width_df['none CTCF'], alternative="greater")[1]

            results.append([celltype, p_value, -np.log10(p_value), withCTCF.shape[0], woCTCF.shape[0]])

        final_df = pd.DataFrame(results)
        final_df.columns = ['celltype', 'p_value', 'negative log10 of p_value', 'number of CTCF TE/SE', 'number of non CTCF TE/SE']
        final_df.to_csv('../CTCF_SE/'+'CTCF_H3K27ac_width_cmp_allcelltypes.xls', sep='\t', index=None)

    """
    Not all SE is associated with CTCF, draw the SE with CTCF height, SE with CTCF distance
    """
    if False:
        celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte', 'cardiac-muscle',
                     'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
                     'neural-cell', 'neutrophil', 'osteoblast', 'smooth-muscle-cell', 'H1-hESC']
        results = []
        for celltype in celltypes:
            CTCF_SE_df = pd.read_csv('../CTCF_SE/'+celltype + '_CTCF_H3K27ac_refmap.xls', sep='\t')

            # cutoff = SE_cutoff(CTCF_SE_df)

            CTCF_SE_df = CTCF_SE_df.nlargest(1000, 'SE_width')
            withCTCF = CTCF_SE_df[CTCF_SE_df['distance'] == 0].copy()
            woCTCF = CTCF_SE_df[CTCF_SE_df['distance'] != 0].copy()

            CTCF_SE_df['distance'] = np.log10(CTCF_SE_df['distance'] + 1)
            CTCF_SE_df['SE_width'] = np.log10(CTCF_SE_df['SE_width'])

            plot = CTCF_SE_df.plot.scatter(x='SE_width', y='distance', edgecolors=None, c='red', s=1)
            fig = plot.get_figure()
            fig.savefig('../CTCF_SE/'+celltype + '_CTCF_H3K27ac_SE_width_vs_distance.pdf')
            plt.close("all")

            width_df = pd.DataFrame(index=range(max(withCTCF.shape[0], woCTCF.shape[0])), columns=['CTCF', 'none CTCF'])
            width_df['CTCF'] = list(withCTCF['SE_width'].values) + [np.nan] * (width_df.shape[0] - withCTCF.shape[0])
            width_df['none CTCF'] = list(woCTCF['SE_width'].values) + [np.nan] * (width_df.shape[0] - woCTCF.shape[0])

            width_df.to_csv('../CTCF_SE/' + celltype + '_CTCF_H3K27ac_SE_width_cmp.xls', sep='\t', index=None)

            p_value = stats.mannwhitneyu(width_df['CTCF'], width_df['none CTCF'], alternative="greater")[1]
            print celltype, width_df['CTCF'].median(), width_df['none CTCF'].median()

            results.append([celltype, p_value, -np.log10(p_value), withCTCF.shape[0], woCTCF.shape[0]])

        final_df = pd.DataFrame(results)
        final_df.columns = ['celltype', 'p_value', 'negative log10 of p_value', 'number of CTCF TE/SE',
                            'number of non CTCF TE/SE']
        final_df.to_csv('../CTCF_SE/' + 'CTCF_H3K27ac_SE_width_cmp_allcelltypes.xls', sep='\t', index=None)

    """
        Not all SE is associated with CTCF, draw the SE with CTCF height, SE with CTCF distance
    """
    if False:
        celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte', 'cardiac-muscle',
                     'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
                     'neural-cell', 'neutrophil', 'osteoblast', 'smooth-muscle-cell']
        for celltype in celltypes:
            CTCF_SE_df = pd.read_csv('../CTCF_SE/' + celltype + '_CTCF_H3K27ac_refmap.xls', sep='\t')

            CTCF_SE_df = CTCF_SE_df[CTCF_SE_df.index.isin(CTCF_SE_df.nlargest(1000, 'SE_width').index)]

            CTCF_SE_df['SE_width'] = np.log10(CTCF_SE_df['SE_width'])
            print celltype

            plot = plt.scatter(x=CTCF_SE_df['SE_width'], y=CTCF_SE_df['CTCF_height'], edgecolors=None, c='red', s=1)
            fig = plot.get_figure()
            fig.savefig('../CTCF_SE/' + celltype + '_CTCF_H3K27ac_SE_width_vs_height.pdf')
            plt.close("all")

    """
        What is the accumulated distrabution of CTCF_SE.
    """
    if False:
        celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte', 'cardiac-muscle',
                     'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
                     'neural-cell', 'neutrophil', 'osteoblast', 'smooth-muscle-cell', 'H1-hESC']
        distances = [1] + range(1000, 101000, 1000)
        result_df = pd.DataFrame(index=distances, columns=celltypes)
        for celltype in celltypes:
            CTCF_SE_df = pd.read_csv('../CTCF_SE/' + celltype + '_CTCF_H3K27ac_refmap.xls', sep='\t')
            CTCF_SE_df['distance'] = np.log10(CTCF_SE_df['distance'] + 1)
            i = 0
            for distance in list(np.log10(distances)):
                cur_df = CTCF_SE_df[CTCF_SE_df['distance']<=distance]
                result_df.ix[distances[i], celltype] = cur_df.shape[0]*1./CTCF_SE_df.shape[0]
                i+=1
        result_df['avg'] = result_df.mean(axis=1)
        result_df['delta_percentage'] = result_df['avg'].diff().shift(-1)
        result_df['distance'] = result_df.index
        result_df['delta_distance'] = result_df.distance.diff().shift(-1)
        result_df['delta_percentage'] = result_df['delta_percentage']/result_df['delta_distance']
        result_df.to_csv('../CTCF_SE/alltypes_cumulative_distance_CTCF_SE.xls', sep='\t')

    """
        CTCF-SE is enriched in SE
    """
    if False:
        SE_df = pd.read_csv('../CTCF_SE/' + 'CTCF_H3K27ac_SE_width_cmp_allcelltypes.xls', sep='\t', index_col=0)
        df = pd.read_csv('../CTCF_SE/' + 'CTCF_H3K27ac_width_cmp_allcelltypes.xls', sep='\t', index_col=0)
        results = []
        for celltype in SE_df.index:
            total = df.ix[celltype, 'number of CTCF TE/SE'] + df.ix[celltype, 'number of non CTCF TE/SE']
            CTCF_SE = SE_df.ix[celltype, 'number of CTCF TE/SE']
            number_SE = 1000
            CTCF_TESE = df.ix[celltype, 'number of CTCF TE/SE']
            p_value = stats.fisher_exact([[CTCF_SE, number_SE-CTCF_SE], [CTCF_TESE-CTCF_SE, total-CTCF_TESE-number_SE+CTCF_SE]])[1]
            results.append((celltype, p_value))
        final_df = pd.DataFrame(results)
        final_df.columns = ['celltype', 'p_value']
        final_df.to_csv('../CTCF_SE/CTCF_SE_enrich_in_SE.xls', sep='\t',index=None)

    """
    CTCF-SE is specific across celltype, SE is specific, but the case is some CTCF-SE become OSE in some celltypes
    """
    # Step 1: merge the super enhancer to get SE reference map
    if False:
        celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte', 'cardiac-muscle',
                         'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
                         'neural-cell', 'neutrophil', 'osteoblast', 'smooth-muscle-cell', 'H1-hESC']
        refmap = None
        for celltype in celltypes:
            cur_df = pd.read_csv('../CTCF_SE/'+celltype + '_CTCF_H3K27ac_refmap.xls', sep='\t')
            cur_df = cur_df.nlargest(3000, 'SE_width').copy()
            if refmap is None:
                refmap = cur_df
            else:
                refmap = refmap.append(cur_df)
        refmap = refmap.sort_values(['SE_chr', 'SE_start'], ascending=[True, True])
        refmap.index = range(refmap.shape[0])
        results = []
        chr, start, end = None, None, None
        for i in range(refmap.shape[0]):
            cur_chr, cur_start, cur_end = refmap.iloc[i, [0,1,2]]
            if cur_chr is None:
                chr, start, end = cur_chr, cur_start, cur_end
            elif cur_chr != chr:
                results.append((chr, start, end))
                chr, start, end = cur_chr, cur_start, cur_end
            elif cur_start - end > 12500:
                results.append((chr, start, end))
                chr, start, end = cur_chr, cur_start, cur_end
            else:
                if cur_end > end:
                    end = cur_end
        results.append((chr, start, end))
        result_df = pd.DataFrame(results)
        result_df.columns =['chr', 'start', 'end']
        result_df.to_csv('../CTCF_SE/H3K27ac_SE_refmap.xls', sep='\t', index=None)

    # Step 2: Call peak in SE and CTCF across 15 different types
    if False:
        celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte', 'cardiac-muscle',
                     'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
                     'neural-cell', 'neutrophil', 'osteoblast', 'smooth-muscle-cell', 'H1-hESC']
        # celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell']
        # celltypes = ['NHLF', 'astrocyte', 'cardiac-muscle',
        #              'monocyte', 'hepatocyte']
        # celltypes = ['fibroblast-dermis', 'myotube', 'keratinocyte',
        #              'neural-cell', 'neutrophil', 'osteoblast', 'smooth-muscle-cell']
        path = '/archive2/tmhbxx3/CancerLearn/data/normal/wigs_qn/'

        #Merge reference map
        refmap = pd.read_csv('H3K27ac_SE_refmap.xls', sep='\t')
        from Wig import *

        for celltype in celltypes:
            cur_wig_path = [x for x in os.listdir(path+celltype+'/H3K27ac/') if x.endswith('.wig')][0]
            cur_wig = Wig(path+celltype+'/H3K27ac/'+cur_wig_path)
            peaks = []
            for i in range(refmap.shape[0]):
                cur_chr, cur_start, cur_end = refmap.ix[i, ['chr', 'start', 'end']]
                if cur_chr not in cur_wig.genome.keys():
                    continue
                cur_peaks = cur_wig.genome[cur_chr].get_peaks(cur_start, cur_end, cutoff=10, min_width=40)
                peaks += cur_peaks

            df = pd.DataFrame(peaks)
            df.columns = ['chr', 'start', 'end', 'center', 'width_above_cutoff', 'total_signal', 'height',
                          'height_logP',
                          'skewness', 'kurtosis']

            df.to_csv(celltype + '_SE_' + str(10) + '.csv', index=None)

            cur_wig_path = [x for x in os.listdir(path + celltype + '/CTCF/') if x.endswith('.wig')][0]
            cur_wig = Wig(path + celltype + '/CTCF/' + cur_wig_path)
            peaks = []
            for i in range(refmap.shape[0]):
                cur_chr, cur_start, cur_end = refmap.ix[i, ['chr', 'start', 'end']]
                if cur_chr not in cur_wig.genome.keys():
                    continue
                cur_peaks = cur_wig.genome[cur_chr].get_peaks(cur_start, cur_end, cutoff=10, min_width=40)
                peaks += cur_peaks

            df = pd.DataFrame(peaks)
            df.columns = ['chr', 'start', 'end', 'center', 'width_above_cutoff', 'total_signal', 'height',
                          'height_logP',
                          'skewness', 'kurtosis']

            df.to_csv(celltype + '_CTCF_' + str(10) + '.csv', index=None)

    # Step 3: to check which one is CSE and which one is OSE and get the specificity
    if False:
        celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte', 'cardiac-muscle',
                     'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
                     'neural-cell', 'neutrophil', 'osteoblast', 'smooth-muscle-cell', 'H1-hESC']
        refmap = pd.read_csv('../CTCF_SE/H3K27ac_SE_refmap.xls', sep='\t')

        SE_df = pd.DataFrame(columns=celltypes)
        CTCF_df = pd.DataFrame(columns=celltypes)

        for celltype in celltypes:
            cur_df = pd.read_csv('../CTCF_SE/specificity/'+celltype+'_SE_10.csv')
            for i in range(refmap.shape[0]):
                cur_chr, cur_start, cur_end = refmap.iloc[i, :3]
                cur_index = cur_chr+':'+str(cur_start)+'-'+str(cur_end)
                cur_SE_df = cur_df[(cur_df.chr==cur_chr)&(cur_df.start>=cur_start)&(cur_df.end<=cur_end)]
                if cur_SE_df.shape[0]==0:
                    cur_width = 0
                    cur_height = 0
                    cur_total_signal = 0
                else:
                    cur_width = cur_SE_df['width_above_cutoff'].sum()
                    cur_total_signal = cur_SE_df['total_signal'].sum()
                    cur_height = cur_SE_df['height'].max()
                SE_df.ix[cur_index, celltype] = cur_width

            cur_df = pd.read_csv('../CTCF_SE/specificity/' + celltype + '_CTCF_10.csv')
            for i in range(refmap.shape[0]):
                cur_chr, cur_start, cur_end = refmap.iloc[i, :3]
                cur_index = cur_chr + ':' + str(cur_start) + '-' + str(cur_end)
                cur_SE_df = cur_df[(cur_df.chr == cur_chr) & (cur_df.start >= cur_start) & (cur_df.end <= cur_end)]
                if cur_SE_df.shape[0] == 0:
                    cur_width = 0
                    cur_height = 0
                    cur_total_signal = 0
                else:
                    cur_width = cur_SE_df['width_above_cutoff'].sum()
                    cur_total_signal = cur_SE_df['total_signal'].sum()
                    cur_height = cur_SE_df['height'].max()
                CTCF_df.ix[cur_index, celltype] = cur_height
        SE_df.to_csv('../CTCF_SE/SE_specificity.xls', sep='\t')
        CTCF_df.to_csv('../CTCF_SE/CTCF_specificity.xls', sep='\t')

    # Check whether CTCF asscociated SE is more specific
    if False:
        SE_df = pd.read_csv('../CTCF_SE/SE_specificity.xls', index_col=0, sep='\t')
        CTCF_df = pd.read_csv('../CTCF_SE/CTCF_specificity.xls', index_col=0, sep='\t')
        SE_max = SE_df.max(axis=1)
        SE_df['tau'] = tau(SE_df)
        CTCF_df['max'] = CTCF_df.max(axis=1)

        CSE_df = SE_df[CTCF_df['max']>=2000]
        OSE_df = SE_df[CTCF_df['max']<50]

        SE_df['max'] = SE_max
        SE_df['CTCF height'] = CTCF_df['max']

        print SE_df[['CTCF height', 'tau']].corr(method='spearman')

        print stats.mannwhitneyu(CSE_df['tau'], OSE_df['tau'], alternative='less')[1]
        print CSE_df['tau'].median(), OSE_df['tau'].median()

        print SE_df[(SE_df['CTCF height']<50) & (SE_df['max']<10000)]

    """
    What is the genes associated with CSE and OSE, we need to remove the genes with CTCF around, see next step
    """
    if False:
        # celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte',
        #              'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
        #              'neural-cell', 'osteoblast', 'smooth-muscle-cell']
        celltypes = ['cardiac-muscle','neutrophil',  'H1-hESC']
        results = []
        for celltype in celltypes:
            CTCF_SE_df = pd.read_csv('../CTCF_SE/' + celltype + '_CTCF_H3K27ac_refmap.xls', sep='\t')

            # cutoff = SE_cutoff(CTCF_SE_df)
            withCTCF = CTCF_SE_df[CTCF_SE_df['distance'] == 0].copy()
            woCTCF = CTCF_SE_df[CTCF_SE_df['distance'] > 10000].copy()

            withCTCF = withCTCF.nlargest(500, 'SE_width')
            woCTCF = woCTCF.nlargest(1500, 'SE_width')

            withCTCF.to_csv('../CTCF_SE/'+celltype+'_CSE.bed', sep='\t', header=None, index=None)
            woCTCF.to_csv('../CTCF_SE/' + celltype + '_OSE.bed', sep='\t', header=None, index=None)

            os.system('/Users/boxia/tools/bedtools2/bin/bedtools intersect -wa -a ../ref_data/hg19.GREATgene2UCSCknownGenes.bed -b '+'../CTCF_SE/'+celltype+'_CSE.bed > '+'../CTCF_SE/'+celltype+'_CSE_genes.xls')
            os.system('/Users/boxia/tools/bedtools2/bin/bedtools intersect -wa -a ../ref_data/hg19.GREATgene2UCSCknownGenes.bed -b ' + '../CTCF_SE/'+celltype + '_OSE.bed > ' + '../CTCF_SE/'+celltype + '_OSE_genes.xls')

            withCTCF = pd.read_csv('../CTCF_SE/'+celltype + '_CSE_genes.xls', sep='\t', header=None)
            woCTCF = pd.read_csv('../CTCF_SE/'+celltype + '_OSE_genes.xls', sep='\t', header=None)
            withCTCF.columns = ['SE_chr', 'SE_start', 'SE_end', 'gene_id']
            woCTCF.columns = ['SE_chr', 'SE_start', 'SE_end', 'gene_id']

            withCTCF_genes = pd.read_excel('../CTCF_SE/'+celltype + '_CSE_genes.xlsx', index_col=0)
            woCTCF_genes = pd.read_excel('../CTCF_SE/' + celltype + '_OSE_genes.xlsx', index_col=0)

            woCTCF = woCTCF[~woCTCF['gene_id'].isin(withCTCF_genes.index)]
            woCTCF = woCTCF[~woCTCF['gene_id'].isin(withCTCF['gene_id'].unique())]

            print len(withCTCF['gene_id'].unique()), len(woCTCF['gene_id'].unique())
            withCTCF.to_csv('../CTCF_SE/'+celltype + '_CSE_genes.xls', sep='\t', index=None)
            woCTCF.to_csv('../CTCF_SE/' + celltype + '_OSE_genes.xls', sep='\t', index=None)


    """
       Get the CTCF associated genes and other genes
    """
    if False:
        # celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte',
        #              'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
        #              'neural-cell', 'osteoblast', 'smooth-muscle-cell']
        celltypes = ['cardiac-muscle','neutrophil', 'H1-hESC']
        results = []
        for celltype in celltypes:
            celltype_df = pd.read_csv('../real_tables/' + celltype + '_parameters_default_real_table.csv', index_col=0)
            celltype_df['CTCF_max'] = celltype_df[['CTCF_height_genebody', 'CTCF_height']].max(axis=1)

            # celltype_df = celltype_df.nlargest(3000, 'H3K27ac_total_width')
            # print celltype_df[['CTCF_height', 'H3K27ac_total_width']].corr(method='spearman')

            CSE_df = celltype_df[(celltype_df['CTCF_height_genebody']>=50)|(celltype_df['CTCF_height']>=50)]
            OSE_df = celltype_df[(celltype_df['CTCF_height_genebody']<50)&(celltype_df['CTCF_height']<50)]
            # CSE_df = CSE_df.nlargest(500, 'H3K27ac_total_width')
            # OSE_df = OSE_df.nlargest(500, 'H3K27ac_total_width')

            print celltype, CSE_df['H3K27ac_total_width'].shape[0], OSE_df['H3K27ac_total_width'].shape[0]

            CSE_df[['CTCF_max', 'H3K27ac_total_width']].to_excel('../CTCF_SE/'+celltype+'_CSE_genes.xlsx')
            OSE_df[['CTCF_max', 'H3K27ac_total_width']].to_excel('../CTCF_SE/' + celltype + '_OSE_genes.xlsx')

    """
        CTCF have low correlation with H3K27ac in coding genes
    """
    if False:
        celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte',
                     'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
                     'neural-cell', 'osteoblast', 'smooth-muscle-cell']
        # 'cardiac-muscle','neutrophil',
        results = []
        for celltype in celltypes:
            celltype_df = pd.read_csv('../real_tables/' + celltype + '_parameters_default_real_table.csv', index_col=0)

            withCTCF = pd.read_csv('../CTCF_SE/' + celltype + '_CSE_genes.xls', sep='\t')
            woCTCF = pd.read_csv('../CTCF_SE/' + celltype + '_OSE_genes.xls', sep='\t')

            withCTCF =  celltype_df[celltype_df.index.isin(withCTCF['gene_id'].unique())]
            woCTCF = celltype_df[celltype_df.index.isin(woCTCF['gene_id'].unique())]

            final_df = pd.DataFrame(columns=['CSE', 'OSE'], index=range(len(list(withCTCF.index)+list(woCTCF.index))))
            final_df['CSE'] = list(withCTCF['H3K27ac_total_width'].values) + [np.nan]*(final_df.shape[0]-withCTCF.shape[0])
            final_df['OSE'] = list(woCTCF['H3K27ac_total_width'].values) + [np.nan] * (
            final_df.shape[0] - woCTCF.shape[0])
            final_df.to_csv('../CTCF_SE/'+celltype+'_CSE_vs_OSE_genes_SEwidth.xls', index=None, sep='\t')


