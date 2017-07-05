import pandas as pd


analysis = pd.read_pickle('/Users/annabroido/Dropbox/Research/LRTAnalysis/LRTAnalysis/analysis/analysis.p')
rows = analysis[analysis.index.str.contains("Route_Views*")]

## SNAP datasets
fp_gml = "/Volumes/Durodon/gmls/Technological/Communication/nn/Route_Views_AS_graphs_1997-2000_Union_Technological_Communication_nn.gml"
for ind,row in rows.iterrows():
    analysis.loc[ind].fp_gml = fp_gml

## Congress
# Senate
rows = analysis[analysis.index.str.contains("Senate")]
fp_gml = "/Volumes/Durodon/gmls/Social/Political/nn/US_Congress_roll-call_votes_1789-2009_Rollcall_US_Senate_Union_Social_Political_nn.gml"
for ind,row in rows.iterrows():
    analysis.loc[ind].fp_gml = fp_gml
# House
rows = analysis[analysis.index.str.contains("House")]
fp_gml = "/Volumes/Durodon/gmls/Social/Political/nn/US_Congress_roll-call_votes_1789-2009_Rollcall_US_House_Union_Social_Political_nn.gml"
for ind,row in rows.iterrows():
    analysis.loc[ind].fp_gml = fp_gml


analysis.to_pickle('/Users/annabroido/Dropbox/Research/LRTAnalysis/LRTAnalysis/analysis/analysis.p')
