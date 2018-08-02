import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

df = pd.read_csv('test.txt', header=None)
print df.iloc[:, 0].unique()