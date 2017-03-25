import numpy as np
import fit
import scipy.optimize as op
import scipy.special as sp
import integration_constants as ic
from scipy.stats import norm

def pllogpdf(x,alpha):
    logpdf = np.log(ic.plconst(np.array([alpha]), np.min(x))) -alpha*np.log(x)
    return logpdf

def vuong(LplV, LaltV):
    logratioV = LplV - LaltV
    n = len(logratioV)
    R = np.sum(logratioV)
    # mean of normal dist
    mu = np.mean(logratioV)
    # standard deviation of normal dist
    sigma = np.std(logratioV)
    normR = (1/np.sqrt(n))*R/sigma
    vuongstat = np.sqrt(n)*mu/sigma
    # one-sided p-value
    p1 = norm.cdf(vuongstat)
    if p1 > 0.5:
        p1 = 1-p1
    # 2-sided p-value
    p2 = 2*p1
    return R, p2, normR


def decide(logratio, p, decisionthresh=0.1):
    if p <= decisionthresh:
        if logratio == 0:
            d = 0
        elif logratio > 0:
            # PL better
            d = -1
        else:
            # exp better
            d = 1
    else:
        # inconclusive
        d = 0
    return d


"""
Perform likelihood ratio test for exponetial distribution. First fits an
exponential distribution to the data. A Vuong statistic is calculated from
the likelihoods of the exponential fit and the power law fit. The
convergence criteria for exponential fitting is met by either hitting the
maximum number of iterations, or by reaching a likelihood value at which a
decision is statistically possible.

Input:
    x           ndarray, data set to be fit
    alpha       float, exponent of power law fit

Output:
    dexp    int, decision about exponential distribution
"""

# define log pdf, for use in likelihood calculation
def logpdf(x,lam):
    xmin = np.min(x)
    ntail = len(x)
    result = np.log(1-np.exp(-lam))+lam*xmin - lam*x
    return result




terrfp = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/terrorism.txt'
wordsfp = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/uniquewords.txt'
xfull = np.loadtxt(terrfp, dtype=int)

[alpha,xmin, ntail, Lpl, ks] = fit.pl(xfull)
x = xfull[xfull>=xmin]
ntail = len(x)
lam0 = np.log(1+float(ntail)/np.sum(x-xmin))
# define negative log likelihood, the function we wish to minimize
negloglike = lambda lam: -np.sum(logpdf(x,lam))
tol = 1E-9
res = op.minimize(negloglike,lam0, bounds=[(tol,None)], method='L-BFGS-B')
lam = np.asscalar(res.x)
loglike = -negloglike(lam)
# perform lrt: Log-likelihood ratio between discrete power law and
# exponential distribution. This is done pointwise so that we can use
# Vuong's statistic to estimate the variance in the ratio
LplV = pllogpdf(x,alpha)
LexpV = logpdf(x,lam)
R, p, normR = vuong(LplV, LexpV)
# check if statistically significant
dexp = decide(R, p)
dexp, R, p, loglike, normR
