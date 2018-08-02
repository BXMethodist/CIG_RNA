import os

genebodys = ['TSS', 'TTS']
markers = ['h3k4me3_qn', 'h3k4me1_qn', 'h3k27me3_qn', 'h3k27ac_qn']
markers = ['h3k4me1_qn']
features = [ 'single_width', 'total_width', 'total_signal', 'single_signal', 'height','skewness', 'kurtosis'  ]
for TTS in genebodys:
    for marker in markers:
        for feature in features:
            f = open(marker +feature+TTS + '.pbs', 'w')
            f.write('#!/bin/bash\n#PBS -r n\n#PBS -m e\n#PBS -M bxia@houstonmethodist.org\n')
            f.write('#PBS -l walltime=96:00:00\n')
            f.write('#PBS -l nodes=1:ppn=8\n')
            f.write('##PBS -q mediummem\n')
            f.write('#PBS -q default\n')
            f.write('##PBS -l pmem=16000mb\n')
            f.write('module load python/2.7.11\n')
            f.write('module load R/3.2.1\n')
            f.write('cd ' + os.getcwd() + '\n')
            if feature != 'skewness' and feature != 'kurtosis':
                cmd = 'python /home/tmhbxx3/archive/tools/danpos_dev/danpos.py grid --TTS_pos ' + TTS + ' --up_stream_grid start:-250000:1000:1000:25000:2:1000 --down_stream_grid end:0:251000:1000:25000:2:1000 --height_grid height:1:41:1:4:2:1 -g 11 -n gene -t cell_type CIG_Bo_curated.xlsx non_CIG_control.csv /home/tmhbxx3/scratch/CIG/' + marker + '_peaks/pooled/ ' + feature + ' ' + marker + '_' + TTS + ' ./ /archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls'
            else:
                cmd = 'python /home/tmhbxx3/archive/tools/danpos_dev/danpos.py grid --TTS_pos ' + TTS + ' --up_stream_grid start:-250000:1000:1000:25000:2:1000 --down_stream_grid end:0:251000:1000:25000:2:1000 --height_grid height:1:41:1:4:2:1 -g 11 -n gene -t cell_type CIG_Bo_curated.xlsx non_CIG_control.csv /home/tmhbxx3/scratch/CIG/test/' + marker + '_sk_peaks/ ' + feature + ' ' + marker + '_' + TTS + ' ./ /archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls'
            f.write(cmd+'\n')
            f.close()
            os.system('qsub ' + marker +feature+TTS + '.pbs')
        #     break
        # break


