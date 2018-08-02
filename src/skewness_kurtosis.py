from Wig import Wig
import os, sys
import numpy as np




if __name__ == "__main__":
    index = int(sys.argv[1])

    path = '/home/tmhbxx3/archive/CIG_fig7/'
    wigs = [x for x in os.listdir(path) if x.endswith('.wig')]

    cur_wig = Wig(path+wigs[index])

    cur_name = wigs[index][:-4]
    cur_out = '/archive/tmhbxx3/CIG_fig7_sk/'+cur_name

    # for cutoff in range(1, 100, 1):
    for cutoff in [31,6,96,48,5,24,13]:
        cur_wig.CallPeak(cutoff, output=cur_out)



