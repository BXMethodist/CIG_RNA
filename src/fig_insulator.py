import pandas as pd, numpy as np, os

"""
whether CTCF SE are close to insulator?
"""

if __name__ == "__main__":
    """
    check human
    """
    if False:
        insulator_df = pd.read_excel('../insulator/human_stem_cell.xlsx', header=1)
        results = []
        for i in range(insulator_df.shape[0]):
            chr1, start1, end1, chr2, start2, end2 = insulator_df.ix[i, :]
            "Chromosome      Start        End Chromosome.1    Start.1      End.1"
            results.append((chr1, start1, end1))
            results.append((chr2, start2, end2))
        result_df = pd.DataFrame(results)
        result_df.columns = ['chr', 'start', 'end']
        result_df = result_df.sort_values(by=['chr', 'start'])
        result_df = result_df.drop_duplicates()
        result_df.to_excel('../insulator/human_stem_cell_unique.xlsx', index=None)

    """
    How many CTCF-SE (SE is defined by top 1000 SE, with CTCF) are in the insulator
    """
    if True:
        pass



