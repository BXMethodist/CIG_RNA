import os, sys, pandas as pd

sys.path.append('/archive/tmhbxx3/tools/tools/')
from utils import *
from selector_utils import *

def RunChipSeq(geo_sample_id):
    cmds = fastq_dump(geo_sample_id)
    cmds += bowtie(geo_sample_id)
    submit_pbs(cmds, geo_sample_id[0])
    return

if __name__ == "__main__":
    if False:
        df = pd.read_csv('data.txt', sep='\t', index_col=0)
        for index in df.index:
            RunChipSeq(index)

    if False:
        cmds = danpos_input('SRR385628', 'SRR385630')
        submit_pbs(cmds, 'SRR385628')

        cmds = danpos_input('SRR385631', 'SRR385633')
        submit_pbs(cmds, 'SRR385631')

        cmds = danpos_input('SRR385634', 'SRR385637')
        submit_pbs(cmds, 'SRR385634')

        cmds = danpos_input('SRR385638', 'SRR385641')
        submit_pbs(cmds, 'SRR385638')

    if False:
        # wigs = [x for x in os.listdir('.') if x.endswith('.wig')]
        # wigs = ['osteoblast_H3K27me3_samples.bgsub.Fnor.wig']
        wigs = ['']
        for wig in wigs:
            cmds = danpos_wig_qn(wig, 'H3K4me3')
            submit_pbs(cmds, wig[:-4], mem="50000mb", queue='highmem', ppn=1, walltime="32:00:00")

    if False:
        cutoffs = [5]
        markers = [x for x in os.listdir('.') if os.path.isdir(x)]

        for marker in markers:
            celltypes = [marker + '/' + x + '/' for x in os.listdir(marker) if os.path.isdir(marker + '/'+x)]
            for celltype in celltypes:
                print celltype
                wigs = [celltype+x for x in os.listdir(celltype) if x.endswith('.wig')]
                print wigs
                for w in wigs:
                    cmds = danpos_peaks(w, cutoffs)
                    submit_pbs(cmds, w.replace('/','_'), mem="16000mb", ppn=1)

    if False:
        # markers = [x for x in os.listdir('.') if os.path.isdir(x)]

        celltypes = ['H1-hESC']

        for celltype in celltypes:
            markers = [celltype+'/'+ x + '/' for x in os.listdir(celltype)]
            for marker in markers:
                cmd = danpos_sk(0, celltype+'/'+marker, os.getcwd() + '/')
                # print cmd
                submit_pbs(cmd, celltype+'/'+marker, mem="16000mb", ppn=1)

    if False:
        wigs = [x for x in os.listdir('.') if x.endswith('.wig')]
        for i in range(len(wigs)):
            cmds = danpos_sk(i, '/home/tmhbxx3/archive2/CancerLearn/data/melonoma/wigs_qn/', '/home/tmhbxx3/archive2/CancerLearn/data/melonoma/wigs_qn/')
            submit_pbs(cmds, wigs[i][:-4]+'_sk', mem="16000mb", queue='mediummem', ppn=1, walltime="32:00:00")

    if False:
        celltypes = [x for x in os.listdir('.') if os.path.isdir(x)]
        print celltypes
        for celltype in celltypes:
            markers = [y for y in os.listdir(celltype) if os.path.isdir(celltype+'/'+y)]
            print markers
            for marker in markers:
                wigs = [z for z in os.listdir(celltype+'/'+marker+'/') if z.endswith('.wig')]
                # print wigs
                for wig in wigs:
                    cmd = wigTobw(celltype+'/'+marker+'/'+wig, 'bws/')
                    # print cmd
                    submit_pbs([cmd], wig[:-4], mem="4000mb", queue='default', ppn=1, walltime="32:00:00")

    if True:
        celltypes = [x for x in os.listdir('.') if os.path.isdir(x)]
        celltypes = ['H1-hESC']
        cutoffs = range(100)
        cutoffs = [5, 10]
        print celltypes
        for celltype in celltypes:
            markers = [y for y in os.listdir(celltype) if os.path.isdir(celltype+'/'+y)]
            print markers
            for marker in markers:
                wigs = [z for z in os.listdir(celltype+'/'+marker+'/') if z.endswith('.wig')]
                # print wigs
                for wig in wigs:
                    cmds = danpos_peaks(celltype+'/'+marker+'/'+wig, cutoffs)
                    # print cmds
                    submit_pbs(cmds, 'peaks_'+wig[:-4], mem="8000mb", queue='mediummem', ppn=1, walltime="96:00:00")

    if True:
        celltypes = [x for x in os.listdir('.') if os.path.isdir(x)]
        celltypes = ['H1-hESC']
        for celltype in celltypes:
            markers = [y for y in os.listdir(celltype) if os.path.isdir(celltype+'/'+y)]
            print markers
            for marker in markers:
                wigs = [z for z in os.listdir(celltype+'/'+marker+'/') if z.endswith('.wig')]
                # print wigs
                for wig in wigs:
                    cmds = danpos_sk(0, celltype+'/'+marker+'/',
                                     'peaks_sk/')
                    # print cmds
                    submit_pbs(cmds, 'sk_'+wig[:-4], mem="8000mb", queue='mediummem', ppn=1, walltime="96:00:00")


    if False:
        # samples, inputs = encode_metadata('hepatocyte_metadata.tsv')
        # samples, inputs = encode_metadata('B_cell_metadata.tsv')
        # samples, inputs = encode_metadata('../osteoblast_metadata.tsv')
        # for target in ['H3K27me3-human']:
        # for target in ['H3K4me3-human']:
        # for target in ['H3K4me1-human']:
        # for target in ['H3K27ac-human']:
        # for target in ['CTCF-human']:
        # for target in ['H3K9me3-human']:
        for target in ['H3K79me2-human', 'H3K9me3-human', 'CTCF-human']:
        # for target in ['H3K79me2-human', 'H3K9me3-human']:
        # for target in ['H3K79me2-human', 'H3K9me3-human', 'CTCF-human',
        #                'H3K4me3-human', 'H3K4me1-human', 'H3K27me3-human', 'H3K27ac-human',]:
        # for target in ['H3K79me2-human', 'H3K9me3-human', 'CTCF-human',
        #                'H3K27me3-human']:
            target_folder = target.replace('-human', '')
            os.system('mkdir '+ target_folder)
            # samples, inputs = encode_metadata('astrocyte_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('neutrophil_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('osteoblast_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('hepatocyte_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('B_cell_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('neural_cell_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('NHLF_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('HSMM_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('HMEC_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('HUVEC_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('fibroblast_dermis_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('keratinocyte_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('foreskin_fibroblast_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('monocyte_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('foreskin_melanocyte_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('cardiac_muscle_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('myotube_metadata.tsv', target=target)
            # samples, inputs = encode_metadata('smooth_muscle_cell_metadata.tsv', target=target)
            samples, inputs = encode_metadata('H1-hESC_metadata.tsv', target=target)


            for sample in samples:
                cmds = []
                cmds += encode_download(sample)
                cmds += bowtie(sample)
                submit_pbs(cmds, target_folder+'/'+sample)
            for input in inputs:
                cmds = []
                cmds += encode_download(input)
                cmds += bowtie(input)
                submit_pbs(cmds, target_folder+'/'+input)

    if False:
        # targets = ['H3K79me2-human', 'H3K9me3-human', 'CTCF-human',
        #                'H3K4me3-human', 'H3K4me1-human', 'H3K27me3-human', 'H3K27ac-human', ]
        targets = ['H3K79me2-human', 'H3K9me3-human', 'CTCF-human']

        # folders = [x for x in os.listdir('.') if os.path.isdir(x) if x.find('CD4')==-1]
        # folders = ['cardiac_muscle', 'keratinocyte', 'monocyte', 'myotube', 'smooth_muscle_cell']
        folders = ['H1-hESC']
        for celltype in folders:
            cur_targets = [y for y in os.listdir(celltype) if os.path.isdir(celltype+'/'+y) and y+'-human' in targets]
            for target in cur_targets:
                target_folder = celltype+'/'+ target
                # samples, inputs = encode_metadata('astrocyte_metadata.tsv', target=target)
                # samples, inputs = encode_metadata('neutrophil_metadata.tsv', target=target)
                # samples, inputs = encode_metadata('osteoblast_metadata.tsv', target=target)
                # samples, inputs = encode_metadata('hepatocyte_metadata.tsv', target=target)
                # samples, inputs = encode_metadata('B_cell_metadata.tsv', target=target)
                # samples, inputs = encode_metadata('neural_cell_metadata.tsv', target=target)
                # samples, inputs = encode_metadata('NHLF_metadata.tsv', target=target)
                # samples, inputs = encode_metadata('HSMM_metadata.tsv', target=target)
                # samples, inputs = encode_metadata('HMEC_metadata.tsv', target=target)
                # samples, inputs = encode_metadata('HUVEC_metadata.tsv', target=target)
                # try:
                samples, inputs = encode_metadata(celltype+'/'+celltype+'_metadata.tsv', target=target+'-human')
                # except:
                #     continue
                os.system('mkdir '+target_folder+'/samples')
                os.system('mkdir '+target_folder+'/inputs')
                print samples
                for sample in samples:
                    cmd = 'mv '+target_folder+'/'+sample+'.bowtie ' + target_folder+'/'+'samples'
                    os.system(cmd)
                for input in inputs:
                    cmd = 'mv '+target_folder+'/'+input+'.bowtie ' + target_folder+'/'+'inputs'
                    os.system(cmd)
                print target_folder
                cmds = danpos_input(target_folder+'/'+'samples', target_folder+'/'+'inputs')
                submit_pbs(cmds, target_folder+'/'+'samples', mem="4000mb", queue='default', ppn=1,)

    if False:
        # targets = ['H3K79me2-human', 'H3K9me3-human', 'CTCF-human',
        #                'H3K4me3-human', 'H3K4me1-human', 'H3K27me3-human', 'H3K27ac-human', ]
        targets = ['H3K79me2-human', 'H3K9me3-human', 'CTCF-human']
        # folders = [x for x in os.listdir('.') if os.path.isdir(x)]
        folders = ['H1-hESC']
        for celltype in folders:
            cur_targets = [y for y in os.listdir(celltype) if os.path.isdir(celltype+'/'+y) and y+'-human' in targets]
            for target in cur_targets:
                target_folder = celltype+'/'+ target
                os.system('mv '+target_folder+'/samples/pooled/*.wig wigs_qn')
        wigs = ['wigs_qn/'+ x for x in os.listdir('wigs_qn/')]
        for w in wigs:
            name = w.replace('archive_tmhbxx3_CancerLearn_data_normal_','')
            os.system('mv '+w + ' '+name)

    if False:
        # targets = ['H3K79me2-human', 'H3K9me3-human', 'CTCF-human',
        #                'H3K4me3-human', 'H3K4me1-human', 'H3K27me3-human', 'H3K27ac-human', ]
        # folders = [x for x in os.listdir('.') if os.path.isdir(x)]
        targets = ['H3K79me2-human', 'H3K9me3-human', 'CTCF-human']
        folders = ['H1-hESC']
        for celltype in folders:
            cur_targets = [y for y in os.listdir(celltype) if os.path.isdir(celltype+'/'+y) and y+'-human' in targets]
            for target in cur_targets:
                target_folder = celltype+'/'+ target
                os.system('rm '+target_folder+'/*.fastq')

    if False:
        targets = ['H3K79me2-human', 'H3K9me3-human', 'CTCF-human',
                       'H3K4me3-human', 'H3K4me1-human', 'H3K27me3-human', 'H3K27ac-human', ]

        files = ['wigs_qn/'+x for x in os.listdir('wigs_qn')]
        for cur_target in targets:
            cur_target = cur_target.replace('-human', '')
            if not os.path.isdir('wigs_qn/'+cur_target):
                os.system('mkdir '+'wigs_qn/'+cur_target)
            for f in files:
                if f.find('_'+cur_target+'_')!=-1:
                    os.system('mv '+f + ' wigs_qn/'+cur_target+'/')

    if False:
        folders = ['CTCF',  'H3K9me3', 'H3K79me2']
        # celltypes = ['HMEC_CTCF_samples.bgsub.Fnor.wig', 'cardiac_muscle_CTCF_samples.bgsub.Fnor.wig', 'CD4_CTCF.wig']
        # celltypes = ['smooth_muscle_cell_CTCF_samples.bgsub.Fnor.wig', 'astrocyte_CTCF_samples.bgsub.Fnor.wig', 'neutrophil_CTCF_samples.bgsub.Fnor.wig']
        # folders = ['CTCF']
        reference = 'B_cell_CTCF_samples.bgsub.Fnor.wig'
        for f in folders:
            wigs = [y for y in os.listdir('wigs/') if y.endswith('.wig') if y.find('H1')!=-1 and y.find(f)!=-1]
            # wigs = [y for y in celltypes if y.endswith('.wig')]
            # print wigs[0]
            for i in range(len(wigs)):
                outprefix = wigs[i].split('_')[0]
                cmds = danpos_wiq(os.getcwd()+'/'+'wigs'+'/'+wigs[i], os.getcwd()+'/'+'wigs'+'/'+reference.replace('CTCF', f), outname=os.getcwd()+'/'+'wigs/'+'wiqs_'+outprefix+'/')
                # cmds = danpos_wiq(os.getcwd() + '/' + cur_f + '/' + wigs[i], os.getcwd() + '/' + cur_f + '/' + reference,
                #                   outname=os.getcwd() + '/' + cur_f + '/wiqs/')
                submit_pbs(cmds, 'wigs'+'/'+wigs[i].replace('_samples.bgsub.Fnor.wig', ''), mem="50000mb", queue='highmem', ppn=1, walltime="32:00:00")


    # Move danpos peaks
    if False:
        # Move sk peaks
        import os

        files = os.listdir('.')
        for f in files:
            info = f.split('_')
            folder = info[0]
            if os.path.isdir(folder):
                pass
            else:
                os.system('mkdir ' + folder)
            marker = info[1]
            if not os.path.isdir(folder + '/' + marker):
                os.system('mkdir ' + folder + '/' + marker)
            name = '_'.join(info[:2] + info[-1:])
            os.system('cp ' + f + ' ' + folder + '/' + marker + '/' + name)
            # break

        # Move peaks
        import os

        files = os.listdir('.')
        for f in files:
            info = f.split('_')
            folder = info[0]
            if os.path.isdir(folder):
                pass
            else:
                os.system('mkdir ' + folder)
            marker = info[1]
            if not os.path.isdir(folder + '/' + marker):
                os.system('mkdir ' + folder + '/' + marker)
            name = '_'.join(info[2:4] + info[-1:])
            os.system('cp ' + f + ' ' + folder + '/' + marker + '/' + name)


    # Print out the tracks
    if False:
        markers = ['H3K4me3', 'H3K4me1', 'H3K27ac', 'H3K27me3', 'CTCF', 'H3K79me2', 'H3K9me3']
        colors = ['0,0,0', '255,0,0','100,100,100','0,0,255','50,150,50','0,255,255','255,0,255',]

        for i in range(len(markers)):
            marker = markers[i]
            color = colors[i]

            get_track_own('CancerLearn/normal/', color, '100', local_path='./', suffix=marker+'.bw')



    # Get the real table
    if False:
        celltypes = os.listdir('/archive2/tmhbxx3/CancerLearn/data/normal/wigs_qn/peaks/')
        for i in range(len(celltypes)):
            cmd = 'python table.py '+str(i)
            submit_pbs([cmd], celltypes[i])


    # selector
    if False:
        selector_df = 'GSM2636047_H3K4me3.regions_gene.xls'
        # gtf_df = '/archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.exon.anno.gtf'
        gtf_df = '/archive/tmhkxc48/ref_data/mm9/mm9.20150218.knownGene.exon.anno.gtf'
        peak_df = 'GSM2636047.H3K4me3.regions_0.0.xls'

        transidToGeneName(selector_df, gtf_df)
        getWidth(selector_df, peak_df)









