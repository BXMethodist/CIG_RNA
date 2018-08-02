import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('../matrix/ten_cell_types_lncRNA.xls', sep='\t', index_col=0)

if False:
    for i in range(df.shape[1]):
        print df.columns[i]
        if df.columns[i].find('width') != -1 or df.columns[i].find('signal') != -1 or df.columns[i].find('height') != -1:
            df[df.columns[i]] = np.log10(df[df.columns[i]]+1)
        cur_hist = df[df.columns[i]].hist(bins=100)
        fig = cur_hist.get_figure()
        fig.savefig('../plots/'+ df.columns[i]+'.pdf')

a = plt.scatter(list(df['h3k4me3_qn_total_width']), list(df['h3k27me3_qn_total_width']))

plt.plot()
plt.show()