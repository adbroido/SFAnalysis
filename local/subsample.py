import numpy as np
import pandas as pd


analysis = pd.read_pickle("/Users/annabroido/Dropbox/Research/LRTAnalysis/LRTAnalysis/analysis/analysis.p")

#subsample by gml path
unique_datasets = np.unique(analysis.fp_gml)
subanalysis = pd.DataFrame(columns = analysis.columns)
for i, dataset in enumerate(unique_datasets):
    query = "fp_gml == '%s'" %dataset
    rows = analysis.query(query)
    num_rows = len(rows)
    ind = int(np.floor(np.random.rand()*num_rows))
    subanalysis.loc[rows.iloc[ind].name] = rows.iloc[ind]


subanalysis.to_csv('/Users/annabroido/Dropbox/Research/LRTAnalysis/LRTAnalysis/analysis/subanalysis.csv')
