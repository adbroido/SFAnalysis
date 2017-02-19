import glob
import numpy as np
import pandas as pd
import fit
import importfiles as imp

"""
Assumes  we've already got the degree sequence files.
"""


# files to import
degseqdirpath = '/Users/annabroido/Dropbox/Research/LRTAnalysis/degreesequences/'
# import the .txt files
degseqfiles = glob.glob(degseqdirpath+'*.txt')
# read in next degree sequence
i = 2 # temporily just pick one of them
degseqfile = degseqfiles[i]
[domain, subdomain, m, gsize, x]  = imp.readdata(degseqfile)
# pull out the name for reference
splitstr = degseqfile.split('/')
name = splitstr[len(splitstr)-1]
name = name[0:len(name)-len('distribution.txt')]
# fit the power law and find xmin, alpha, and the log likelihood of the fit
[alpha, xmin, Lpl] = fit.pl(x)
# truncate the data
x = x[x>=xmin]
# compare the alternative distributions
# compare the non-nested alternatives, return the decisions for each
[dexp, dln, dstrexp] = lrt.nonnested(x, Lpl)
# fit the nested alternatives
[dpls, dplwc] = lrt.nested(x, Lpl)
# save the output
analysis = pd.Series([name, xmin, dexp, dln, dstrexp, dpls, dplwc])
