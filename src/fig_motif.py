import pandas as pd, os, sys, numpy as np

sys.path.append('/archive/tmhbxx3/tools/tools/')
from utils import *

if __name__ == "__main__":
    # Run danpos profile for CSE.bed and OSE.bed
    if False:
        celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte', 'cardiac-muscle',
                     'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
                     'neural-cell', 'neutrophil', 'osteoblast', 'smooth-muscle-cell', 'H1-hESC']

        hg19_sizes = pd.read_csv('/archive/tmhbxx3/ref_data/hg19/hg19_chr_sizes.txt', sep='\t', header=None)
        hg19_sizes.columns = ['chr', 'size']
        chr_index = np.random.choice(hg19_sizes.shape[0], size=500, replace=True)
        random_control = open('control.bed', 'w')
        for random_chr in chr_index:
            cur_chr = hg19_sizes.ix[random_chr, 'chr']
            cur_size = hg19_sizes.ix[random_chr, 'size'] - 10000
            cur_start = np.random.choice(cur_size)
            random_control.write(cur_chr+'\t'+str(cur_start)+'\t'+str(cur_start+10000)+'\n')
        random_control.close()

        for celltype in celltypes:
            cmds = ['python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py profile /archive/tmhbxx3/tools/TFmotif/TFBS.wig --bed3file_paths /archive2/tmhbxx3/CancerLearn/data/normal/motif/'+celltype+'_CSE.bed,/archive2/tmhbxx3/CancerLearn/data/normal/motif/'+celltype+'_OSE.bed,/archive2/tmhbxx3/CancerLearn/data/normal/motif/control.bed --wigfile_aliases TFBS --bed3file_aliases CSE,OSE,control --name /archive2/tmhbxx3/CancerLearn/data/normal/motif/'+celltype.replace('-','').replace('_','')+'TF --genomic_sites center --heatmap 1 --flank_up 10000 --flank_dn 10000 --plot_colors red,blue,black']
            cmds.append('python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py profile /archive/tmhbxx3/tools/TFmotif/motif.wig --bed3file_paths /archive2/tmhbxx3/CancerLearn/data/normal/motif/'+celltype+'_CSE.bed,/archive2/tmhbxx3/CancerLearn/data/normal/motif/'+celltype+'_OSE.bed,/archive2/tmhbxx3/CancerLearn/data/normal/motif/control.bed --wigfile_aliases motif --bed3file_aliases CSE,OSE,control --name /archive2/tmhbxx3/CancerLearn/data/normal/motif/'+celltype.replace('-','').replace('_','') +'motif --genomic_sites center --heatmap 1 --flank_up 10000 --flank_dn 10000 --plot_colors red,blue,black')
            submit_pbs(cmds, celltype, mem="8000mb", queue='default', ppn=1, walltime="32:00:00")

    # CSE and OSE signal on H3K9me3, H3K4me1, H3K27me3
    if False:
        celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte',
                     'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
                     'neural-cell', 'osteoblast', 'smooth-muscle-cell', 'H1-hESC']
        celltypes += ['cardiac-muscle', 'neutrophil']
        # markers = ['H3K27me3', 'H3K4me1', 'H3K4me3', 'H3K79me2', 'H3K9me3']
        # markers = ['H3K9me3']
        markers = ['H3K27ac']
        for celltype in celltypes:
            for marker in markers:
                wig = [x for x in os.listdir('/archive2/tmhbxx3/CancerLearn/data/normal/wigs_qn/'+celltype+'/'+marker+'/') if x.endswith('.wig')][0]
                cmds = [
                'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py profile /archive2/tmhbxx3/CancerLearn/data/normal/wigs_qn/'+celltype+'/'+marker+'/'+wig+' --bed3file_paths /archive2/tmhbxx3/CancerLearn/data/normal/motif/' + celltype + '_CSE.bed,/archive2/tmhbxx3/CancerLearn/data/normal/motif/' + celltype + '_OSE.bed,/archive2/tmhbxx3/CancerLearn/data/normal/motif/control.bed --wigfile_aliases TFBS --bed3file_aliases CSE,OSE,control --name /archive2/tmhbxx3/CancerLearn/data/normal/motif/' + celltype.replace(
                    '-', '').replace('_',
                                     '') + marker+' --genomic_sites center --heatmap 1 --flank_up 10000 --flank_dn 10000 --plot_colors red,blue,black']

                submit_pbs(cmds, celltype+marker, mem="8000mb", queue='default', ppn=1, walltime="32:00:00")

    # CSE genes and OSE genes expression different?


    # CSE genes and OSE genes H3K4me1, H3K27ac, H3K4me3 profile
    if False:
        celltypes = ['CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte',
                     'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
                     'neural-cell', 'osteoblast', 'smooth-muscle-cell', 'H1-hESC']
        celltypes += ['cardiac-muscle', 'neutrophil']
        # markers = ['H3K27me3', 'H3K4me1', 'H3K4me3', 'H3K79me2', 'H3K9me3']
        # markers = ['H3K9me3']
        markers = ['H3K27ac']
        gtf_df = pd.read_csv('/home/tmhbxx3/archive/ref_data/hg19/hg19.GREATgene2UCSCknownGenes.table.xls',
                             sep='\t', index_col=0)

        gtf_df.sample(500, random_state=0).to_csv('control_genes.xls', sep='\t', header=None)

        for celltype in celltypes:
            for marker in markers:
                wig = [x for x in os.listdir(
                    '/archive2/tmhbxx3/CancerLearn/data/normal/wigs_qn/' + celltype + '/' + marker + '/') if
                       x.endswith('.wig')][0]
                cse_df = pd.read_csv(celltype + '_CSE_genes.xls', sep='\t')
                ose_df = pd.read_csv(celltype + '_OSE_genes.xls', sep='\t')
                try:
                    cse_df = gtf_df[gtf_df.index.isin(cse_df['gene_id'].unique())]
                    ose_df = gtf_df[gtf_df.index.isin(ose_df['gene_id'].unique())]
                except:
                    cse_df = gtf_df[gtf_df.index.isin(cse_df.iloc[:, 3].unique())]
                    ose_df = gtf_df[gtf_df.index.isin(ose_df.iloc[:, 3].unique())]

                # cse_df.to_csv(celltype + '_CSE_gene.xls', sep='\t', header=None)
                # ose_df.to_csv(celltype + '_OSE_gene.xls', sep='\t', header=None)

                cmds = [
                    'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py profile /archive2/tmhbxx3/CancerLearn/data/normal/wigs_qn/' + celltype + '/' + marker + '/' + wig + ' --genefile_paths /archive2/tmhbxx3/CancerLearn/data/normal/CSE_genes/' + celltype + '_CSE_gene.xls,/archive2/tmhbxx3/CancerLearn/data/normal/CSE_genes/' + celltype + '_OSE_gene.xls,/archive2/tmhbxx3/CancerLearn/data/normal/CSE_genes/control_genes.xls --wigfile_aliases TFBS --genefile_aliases CSE,OSE,control --name /archive2/tmhbxx3/CancerLearn/data/normal/CSE_genes/' + celltype.replace(
                        '-', '').replace('_',
                                         '') + marker + ' --genomic_sites TSS --heatmap 1 --flank_up 3000 --flank_dn 10000 --plot_colors red,blue,black']
                print celltype, marker
                submit_pbs(cmds, celltype + marker, mem="8000mb", queue='default', ppn=1,
                           walltime="32:00:00")

    # Pathway specific

    # differential expressed towards stem cell?


