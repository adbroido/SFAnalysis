import lrt
import numpy as np
import fit

terrfp = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/terrorism.txt'
wordsfp = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/uniquewords.txt'
x = np.loadtxt(terrfp, dtype=int)

[alpha,xmin, ntail, Lpl, ks] = fit.pl(x)
xtail = x[x>=xmin]
dexp, R, p, Lexp, normR = lrt.exp(xtail,alpha, decisionthresh=0.1)
