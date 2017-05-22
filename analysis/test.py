import numpy as np
import fit
import lrt
import pandas as pd
from scipy.stats import norm

fp = "/Users/annabroido/Dropbox/Research/LRTAnalysis/MatlabandRcode/pli-R-v0.0.3-2007-07-25/datasets/uniquewords.txt"
fp = "/Users/annabroido/Dropbox/Research/LRTAnalysis/MatlabandRcode/pli-R-v0.0.3-2007-07-25/datasets/terrorism.txt"
x = np.loadtxt(fp, dtype=int)
[alpha, xmin, ntail, Lpl, ks] = fit.pl(x)
x = x[x>=xmin]
[lam, LexpV, convstatus] = fit.exp(x)
LplV = lrt.pllogpdf(x,alpha)
LaltV = LexpV
logratioV = LplV - LaltV
n = len(logratioV)
R = np.sum(logratioV)
# standard deviation of normal dist
sigma = np.std(logratioV)
normR = (1/np.sqrt(n))*R/sigma
# one-sided p-value
p1 = norm.cdf(normR)
if p1 > 0.5:
    p1 = 1-p1
# 2-sided p-value
p2 = 2*p1
