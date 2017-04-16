import pandas as pd
import numpy as np
import glob

"""
Add the processed individual output files to the main dataframe
"""

df = pd.read_pickle('/Users/anbr3575/LRTAnalysis/analysis/analysis.p')

outputdp = '/Users/anbr3575/LRTAnalysis/output/'
outputfiles = glob.glob(outputdp+'*.csv')
for fp in outputfiles:
    new = pd.read_csv(fp, index_col=0)
    df.loc[new.index[0],list(new.columns)] = list(new.values[0])

df = to_pickle('/Users/anbr3575/LRTAnalysis/analysis/analysis.p')
