import pandas as pd, numpy as np, os
from Wig import Wig
from collections import defaultdict
from bisect import bisect
from scipy import stats

def scale(y, c=True, sc=True):
    x = y.copy()

    if c:
        x -= x.mean()
    if sc and c:
        x /= x.std()
    elif sc:
        x /= np.sqrt(x.pow(2).sum().div(x.count() - 1))
    return x

def peaks_distances(peaks):
    """
    :return:
    """
    distances = []
    target_peak = peaks[0]
    for i in range(1, len(peaks)):
        cur_peak = peaks[i]
        distances.append(cur_peak[0] - target_peak[1])
    return distances

def genome_size(path="TAIR10_chr_sizes.txt", step=10):
    genome = {}
    genome_size_file = open(path, "r")
    for line in genome_size_file.readlines():
        chr_name, chr_size = line.split("\t")
        vector_size = int(chr_size.rstrip())
        if vector_size % step == 0:
            vector_size /= step
        else:
            vector_size = vector_size / step + 1
        genome[chr_name] = np.zeros(vector_size)
    genome_size_file.close()
    return genome

def gini(list_of_values):
    sorted_list = sorted(list_of_values)
    height, area = 0, 0
    for value in sorted_list:
        height += value
        area += height - value / 2.
    fair_area = height * len(list_of_values) / 2.
    return (fair_area - area) / fair_area

def distribution_region_distance(refmap):
    """
    get the distance distribution between the reference map
    :return:
    """
    input_file = open(refmap, "r")

    distance_distributions = []
    chr_name = None

    region_map = defaultdict(list)

    for line in input_file.readlines()[1:]:
        line = line.rstrip().split(',')
        # print line
        chr_name = line[0]
        start = int(line[1])
        end = int(line[2])
        region_map[chr_name].append((start, end))

    for key in region_map.keys():
        region_map[key] = sorted(region_map[key], key=lambda x:x[0])

    for value in region_map.values():
        indexes = [x[0] for x in value]
        for i in range(len(value)-1):
            # cur_start = value[i][0]
            # cur_end_index = bisect(indexes, cur_start + 100000)
            # # cur_peaks = value[i: cur_end_index]
            cur_peaks = value[i: i+2]
            # print cur_peaks
            distance_distributions += peaks_distances(cur_peaks)

    distance_file = open('./10cell_region_gap_distance.txt', 'w')
    for line in distance_distributions:
        # line = [str(l) for l in line]
        distance_file.write(str(line)+'\n')
    distance_file.close()
    return distance_distributions

def get_distribution(file_path, step=100):
    """
    this is function take a txt file, with single columns and return is distribution based on step
    :param file_path:
    :param step:
    :return:
    """
    df = pd.read_csv(file_path, index_col=None, header=None)
    df.columns = ['distance']
    indexes = df['distance'].unique().tolist()
    # indexes = sorted(indexes, key=lambda x: float(x))
    indexes =[x for x in range(100, 100000, 100)]
    # print len(indexes)
    results = []
    counts = {}
    cum = 0
    for i in indexes:
        cur_df = df[df['distance']<=i]
        counts[i] = cur_df.shape[0]
        cum += counts[i]
        p = cur_df.shape[0]*1.0/df.shape[0]
        results.append([i, p])
    result_df = pd.DataFrame(results)
    result_df.to_csv(file_path[:-4]+'distribution.csv', index=None, header=None)
    return result_df

def load_map(reference_map_path, sep=','):
    """
    load a reference map from csv or txt file and change is to dict of lists of list
    :param reference_map_path:
    :return:
    """
    reference_genome_map = {}
    df = pd.read_csv(reference_map_path, sep=sep)
    for chromosome in df.chr.unique():
        cur_df = df[df.chr==chromosome]
        reference_genome_map[chromosome] = cur_df[['start', 'end']].values
    return reference_genome_map

def merge_map(reference_map, distance):
    """
    merge reference map based on the distance provided
    :param reference_map: list of lists containing start, end from same chromosome
    :param distance: cutoff
    :return: merged list of lists
    """
    new_reference_map = []
    start, end = None, None
    for r in reference_map:
        cur_start, cur_end = r
        # print r
        if start is None:
            start, end = cur_start, cur_end
        elif distance >= cur_start - end:
            end = cur_end
        else:
            new_reference_map.append([start, end])
            start = cur_start
            end = cur_end
    if start is not None:
        new_reference_map.append([start, end])
    # print len(new_reference_map)
    return new_reference_map

def get_refmap(cutoffs, outdirectory, genome_size_file='/home/tmhbxx3/archive/ref_data/hg19/hg19_chr_sizes.txt'):
    ## get 10 cell type reference map
    if not outdirectory.endswith('/'):
        outdirectory+='/'
    if not os.path.isdir(outdirectory):
        os.system('mkdir '+ outdirectory)
    for cutoff in cutoffs:
        refmap = refMap(genome_size_file)
        markers = ['h3k4me3_qn', 'h3k4me1_qn', 'h3k27ac_qn', 'h3k27me3_qn']
        for marker in markers:
            peak_files = [f for f in os.listdir('/home/tmhbxx3/scratch/CIG/'+marker+'_peaks/pooled/')
                          if f.endswith('_'+cutoff+'.xls')]
            for peak_file in peak_files:
                refmap.load_signal('/home/tmhbxx3/scratch/CIG/'+marker+'_peaks/pooled/'+peak_file)
        refmap.saveRefMap(outdirectory+'10celltypes_4markers_'+cutoff+'.csv')

def get_map_signal(maps, marker, wig_path, outdirectory):
    wig_path += '/' if not wig_path.endswith('/') else ''
    outdirectory += '/' if not outdirectory.endswith('/') else ''
    if not os.path.isdir(outdirectory):
        os.system('mkdir ' + outdirectory)
    for map in maps:
        cell_refmap = pd.read_csv(map)
        cutoff = int(map.split('_')[-1][:-4])
        index = []
        for i in range(cell_refmap.shape[0]):
            cur_chr, cur_start, cur_end = cell_refmap.ix[i, 'chr'], cell_refmap.ix[i, 'start'], cell_refmap.ix[i, 'end']
            index.append('_'.join([cur_chr, str(cur_start), str(cur_end)]))

        path = wig_path + marker + '/'
        wigs = [x for x in os.listdir(path) if x.endswith('.wig')]
        results_signal = defaultdict(list)
        results_width = defaultdict(list)
        results_height = defaultdict(list)
        for w in wigs:
            cur_wig = Wig(path + w)
            for i in range(cell_refmap.shape[0]):
                cur_chr, cur_start, cur_end = cell_refmap.ix[i, 'chr'], cell_refmap.ix[i, 'start'], cell_refmap.ix[i, 'end']
                cur_signal = cur_wig.genome[cur_chr].get_signals(cur_start, cur_end)[1, :].copy()
                leftover_signal = cur_signal[cur_signal >= cutoff]
                if len(leftover_signal) == 0:
                    cur_total_signal = 0
                    cur_width = 0
                    cur_height = 0
                else:
                    cur_total_signal = leftover_signal.sum()
                    cur_width = len(leftover_signal)
                    cur_height = leftover_signal.max()
                results_signal[w[:-4]].append(cur_total_signal)
                results_width[w[:-4]].append(cur_width)
                results_height[w[:-4]].append(cur_height)
        cur_signal_df = pd.DataFrame.from_dict(results_signal)
        cur_width_df = pd.DataFrame.from_dict(results_width)
        cur_height_df = pd.DataFrame.from_dict(results_height)
        cur_signal_df.index = index
        cur_width_df.index = index
        cur_height_df.index = index

        cur_signal_df.to_csv(outdirectory+marker+'_' + 'signal' + '_'+str(cutoff) +'.csv')
        cur_width_df.to_csv(outdirectory+marker + '_' + 'width' + '_'+str(cutoff) + '.csv')
        cur_height_df.to_csv(outdirectory+marker + '_' + 'height' + '_'+str(cutoff) + '.csv')


def get_map_sd(map):
    df = map.std(axis=1).to_frame()

    df.to_csv(map[:-4]+'_sd.csv')
    return df




class refMap:
    ### take a directory of wig files, generate reference map
    ### based on number of iterations, generate saturation map for coverage, average peak size and number of average peak number
    ### use plotSaturation.py to make figure
    def __init__(self, genome_size_path="TAIR10_chr_sizes.txt"):
        self.genome = genome_size(genome_size_path)
        self.genome_size_path=genome_size_path

    def load_signal(self, path):
        file = open(path, "rb")
        info = file.readlines()
        for line in info:
            # print line
            info = line.split("\t")
            if info[1] == "start":
                continue
            start = int(info[1])/10
            end = int(info[2])/10
            height = float(info[6])
            # print height
            chrName = info[0]
            if chrName in self.genome:
                self.genome[chrName][start - 1:end] = 1
            else:
                print chrName
    def saveRefMap(self, outputname):
        refmap = {}
        for chr, vector in self.genome.items():
            if vector[0] == 1:
                vector = np.insert(vector, 0, 0)
            if vector[-1] == 1:
                vector = np.append(vector, 0)

            # sign change mark the start of the peak and the end to the peak, the end mark is exclusive
            # and this is the index of vector with step = 10, not real genome position
            signchange = ((np.roll(vector, 1) - vector) != 0).astype(int)
            peaksindex = np.where(signchange == 1)[0]

            rowNumber = peaksindex.shape[0] / 2
            colNumber = 2
            peaksindex = peaksindex.reshape((rowNumber, colNumber))

            refmap[chr] = peaksindex

        output = open(outputname + "_refmap.csv", "w")
        output.write('chr,start,end\n')
        for chr, index in refmap.items():
            for i in range(index.shape[0]):
                if (index[i, 1] - index[i, 0]) * 10 > 1000:
                    output.write(chr+','+str(index[i,0]*10)+','+str(index[i,1]*10)+'\n')
        output.close()

if __name__ == "__main__":

    ## get 10 cell type reference map
    # get_refmap([str(x) for x in range(1, 101)], outdirectory='/home/tmhbxx3/scratch/CIG/refmap')

    ## get 10 cell type map signals
    # markers = ['h3k4me3_qn', 'h3k4me1_qn', 'h3k27ac_qn', 'h3k27me3_qn']
    # maps = ['/home/tmhbxx3/scratch/CIG/refmap/'+x for x in os.listdir('/home/tmhbxx3/scratch/CIG/refmap/')]
    # for marker in markers:
    #     get_map_signal(maps=maps, marker=marker, wig_path='/home/tmhbxx3/scratch/CIG/', outdirectory='/home/tmhbxx3/scratch/CIG/heatmap')

    ## Make heatmap
    # df = pd.read_csv('h3k4me3_qn_heatmap5kb.cdt.txt', sep='\t', index_col=0)
    # final = []
    # for i in range(df.shape[0]):
    #     cur_final = df.ix[i, :]
    #     df.ix[i, :] = scale(cur_final)
    # df.to_csv('h3k4me3_qn_heatmap5kb.cdt.txt', sep='\t')

    ## get the overlap region
    df1 = pd.read_csv('h3k4me3_qn_heatmap5kb.csv', index_col=0)
    df2 = pd.read_csv('h3k4me1_qn_heatmap5kb.csv', index_col=0)
    df3 = pd.read_csv('h3k27ac_qn_heatmap5kb.csv', index_col=0)
    df4 = pd.read_csv('h3k27me3_qn_heatmap5kb.csv', index_col=0)
    dfs = [df1, df2, df3, df4]
    markers = ['h3k4me3', 'h3k4me1', 'h3k27ac', 'h3k27me3']
    candidates = {}
    dfs_candidates = {}
    for i in range(len(dfs)):
        df = dfs[i]
        # df['Coefficient of Variation'] = (-df*np.log2(df)).sum(axis=1)
        # print df.max(axis=1)
        # print (1 - df.T/df.max(axis=1)).T
        # lala
        # df['Coefficient of Variation'] = (1./(df.shape[1]-1)) * np.sum((1 - df.T/df.max(axis=1)).T, axis=1)
        df['Coefficient of Variation'] = df.std(axis=1)/df.mean(axis=1)
        # df['Coefficient of Variation'] = (-p*np.log2(p)).sum(axis=1)
        df = df.sort_values(by=['Coefficient of Variation'], ascending=False)
        df['Rank'] = range(1, df.shape[0]+1)
        df.to_csv(markers[i] + '_variation.csv')
        # df = df[df['Coefficient of Variation']>3].copy()
        candidates[markers[i]] = df

    ref_map = pd.read_csv('10_celltypes_4markers_region_refmap.csv')

    l5 = [v[0:3] for v in ref_map.values]

    for i in range(len(l5)):
        l = '_'.join([str(k) for k in l5[i]])
        l5[i] = l
    l5 = [x.split('_') for x in l5]
    # print cur_l[0]
    l5 = [(x[0], int(x[1]), int(x[2])) for x in l5]
    l5 = pd.DataFrame(l5)

    l5.columns = ['chr', 'start', 'end']

    # print l5[0]

    # cig_df = pd.read_excel('/Users/boxia/PycharmProjects/CIG_logistic_regression/CIG_Bo_curated5.xlsx', index_col=0)
    cig_df = pd.read_csv('/Users/boxia/PycharmProjects/CIG_logistic_regression/non_CIG_control_v5.csv', index_col=0)
    gtf_df = pd.read_csv('/Users/boxia/PycharmProjects/CIG_logistic_regression/hg19.ucscgenes.knowngene.xls', sep='\t')

    cig_df = gtf_df[gtf_df['hg19.kgXref.geneSymbol'].isin(cig_df.index)]
    cig_df = cig_df[['hg19.kgXref.geneSymbol', 'hg19.knownGene.chrom', 'hg19.knownGene.txStart', 'hg19.knownGene.txEnd',]]
    cig_df.columns = ['gene', 'chr', 'start', 'end']
    cig_df = cig_df.set_index(['gene'])
    # print cig_df
    # cig_df.to_csv('CIG_transcript.csv')

    # candidates = {'refmap': l5}
    #
    # from itertools import combinations
    #
    # for i in range(2, 4):
    #     combs = combinations(candidates.keys(), i)
    #     for comb in combs:
    #         print comb
    #         result = set()
    #         for r in candidates[comb[0]]:
    #             foo = True
    #             for c in comb[1:]:
    #                 if r not in candidates[c]:
    #                     foo = False
    #                     continue
    #             if foo:
    #                 result.add(r)
    #         print len(result)
    #
    # l5 = [x for x in l1 if x in l2 and x in l3 and x in l4]

    final_selected = set()
    final = {}
    for key in candidates.keys():
        cur_l = candidates[key].index
        cur_l = [x.split('_') for x in cur_l]
        # print cur_l[0]
        cur_l = [(x[0], int(x[1]), int(x[2])) for x in cur_l]
        cur_df = pd.DataFrame(cur_l)

        cur_df.columns = ['chr', 'start', 'end']
        selected = set()
        ranks = defaultdict(int)
        cvs = defaultdict(float)
        for i in range(cig_df.shape[0]):
            cur_chr, cur_start, cur_end = cig_df.ix[i, :]

            # if cig_df.index[i] == 'MYC':
            #     print cur_chr, cur_start, cur_end

            sub_df = cur_df[(cur_df['chr']==cur_chr) &
                            ((cur_df['start'].between(cur_start, cur_end)) |
                             (cur_df['end'].between(cur_start,cur_end)) |
                             ((cur_df['start']<=cur_start) & (cur_df['end']>=cur_end)))]
            if sub_df.shape[0] > 0:
                cur_regions = [value for value in sub_df.values]
                cur_regions_index = []
                for cur_r in cur_regions:
                    cur_r = '_'.join([str(x) for x in cur_r])
                    cur_regions_index.append(cur_r)
                best_rank = candidates[key].ix[cur_regions_index, 'Rank'].min()
                best_cv = candidates[key].ix[cur_regions_index, 'Coefficient of Variation'].max()
                if ranks[cig_df.index[i]] > best_rank or ranks[cig_df.index[i]] == 0:
                    ranks[cig_df.index[i]] = best_rank
                    cvs[cig_df.index[i]] = best_cv
                selected.add(cig_df.index[i])
        print type(ranks)
        ranks_r = []
        for k in ranks.keys():
            ranks_r.append([k, ranks[k], cvs[k]])
        r_df = pd.DataFrame(ranks_r)
        r_df.columns=['gene_id', 'rank', 'sd']
        r_df.to_csv(key+'_nonCIG_ranks.csv')
        print key, len(selected)
        final[key] = selected
        final_selected = final_selected.union(selected)
    print len(final_selected)

    # bo_df = pd.read_excel('/Users/boxia/PycharmProjects/CIG_logistic_regression/CIG_Bo_curated5.xlsx', index_col=0)
    # genes = set(bo_df.index)
    # print len(genes)
    # for gene in genes:
    #     if gene not in final_selected:
    #         cur_cig_df = cig_df[cig_df.index==gene]
    #         # print cur_cig_df
    #         for i in range(cur_cig_df.shape[0]):
    #             cur_chr, cur_start, cur_end = cur_cig_df.iloc[i, 0:3]
    #
    #             # if cig_df.index[i] == 'MYC':
    #             #     print cur_chr, cur_start, cur_end
    #
    #             sub_df = l5[(l5['chr'] == cur_chr) &
    #                             ((l5['start'].between(cur_start, cur_end)) |
    #                              (l5['end'].between(cur_start, cur_end)) |
    #                              ((l5['start'] <= cur_start) & (l5['end'] >= cur_end)))]
    #             if sub_df.shape[0] > 0:
    #                 sub_df['width'] = sub_df['end']-sub_df['start']
    #                 print gene, sub_df
    cig_number = 254
    non_cig_number = 285
    total_regions = 14926
    markers = ['h3k4me3', 'h3k4me1', 'h3k27ac', 'h3k27me3']
    for marker in markers:
        df = pd.read_csv(marker+'_nonCIG_ranks.csv')
        results = []
        for i in np.arange(500, 10000, step=500):
            size = df[(df['rank']>i-500) & (i>df['rank'])].shape[0]
            results.append([i, size, -np.log10(stats.fisher_exact([[size, 254-size], [500-size, total_regions - non_cig_number - 500 + size]],alternative='greater')[1])])
        pd.DataFrame(results).to_csv(marker+'_nonCIG_rank_dis.csv')




