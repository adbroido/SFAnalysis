import sys
sys.path.insert(0, '/Users/annabroido/Dropbox/Research/LRTAnalysis/LRTAnalysis/analysis/')
import pandas as pd
import numpy as np
import importfiles as im

analysis = pd.read_pickle('/Users/annabroido/Dropbox/Research/LRTAnalysis/LRTAnalysis/analysis/analysis.p')
analysis = analysis.query("Domain != 'Economic'")
analysis = analysis.query("alpha != ''")
dsfp = '/Users/annabroido/Dropbox/Research/LRTAnalysis/degreesequences/' # degree sequences

sub = analysis.query('n==""')

unique_datasets = np.unique(sub.fp_gml)
for dataset in (unique_datasets):
    query = "fp_gml == '%s'" %dataset
    rows = analysis.query(query)
    ns = np.zeros(len(rows))
    name = rows.iloc[0].name
    fp = dsfp + name
    x = im.readdata(fp)
    n = len(x)
    ns = n*np.ones(len(rows))
    analysis.loc[analysis.query(query).index,'n']=n
#
analysis.to_pickle('/Users/annabroido/Dropbox/Research/LRTAnalysis/LRTAnalysis/analysis/analysis.p')
analysis.to_csv('/Users/annabroido/Dropbox/Research/LRTAnalysis/LRTAnalysis/analysis/analysis.csv')
