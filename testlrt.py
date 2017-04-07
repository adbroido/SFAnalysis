import numpy as np
import fit
import scipy.optimize as op
import scipy.special as sp
import integration_constants as ic
from scipy.stats import norm, chi2
import lrt

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
            # alt better
            d = 1
    else:
        # inconclusive
        d = 0
    return d

def decidenested(logratio, p, decisionthresh):
    """
    Decisions are   1  -   power law better
                    0  -   inconclusive
                   -1  -   alternative dist better
    """
    if p <= decisionthresh and logratio < 0:
            # alt better
            d = -1
    else:
        # inconclusive
        d = 0
    return d


def ln(x,alpha, decisionthresh):
    """
    Perform likelihood ratio test for log normal distribution. First
    fits a log normal distribution to the data. A Vuong statistic is
    calculated from the likelihoods of the exponential fit and the power law
    fit. The convergence criteria for exponential fitting is met by either
    hitting the maximum number of iterations, or by reaching a likelihood value
    at which a decision is statistically possible.

    Input:
        x           ndarray, data set to be fit
        alpha       float, exponent of power law fit

    Output:
        dstrexp    int, decision about exponential distribution
    """
    def logpdf(x, mu, sigma):
        xmin = np.min(x)
        F = lambda x: (sp.erfc((np.log(x)-mu)/(np.sqrt(2)*sigma)))/2
        g = lambda x: F(x)- F(x+1)
        h = -np.log(F(xmin))+np.log(g(x))
        return h
    # initial estimates
    mu0 = 0
    sigma0 = 1
    theta0 = np.array([mu0, sigma0])
    n = len(x)
    # optimize
    negloglike = lambda theta: -np.sum(logpdf(x,theta[0],theta[1]))
    tol = 1E-1
    bnds=[(-n/5,None),(tol,None)]
    res = op.minimize(negloglike, theta0, bounds=bnds, method='L_BFGS_B')
    theta = np.res.x
    loglike = -negloglike(lam)
    # perform lrt: Log-likelihood ratio between discrete power law and
    # exponential distribution. This is done pointwise so that we can use
    # Vuong's statistic to estimate the variance in the ratio
    LplV = pllogpdf(x,alpha)
    LlnV = logpdf(x,theta[0], theta[1])
    R, p, normR = vuong(LplV, LlnV)
    # check if statistically significant
    dexp = decide(R, p, decisionthresh)
    return dexp, R, p, loglike, normR

terrfp = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/terrorism.txt'
wordsfp = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/uniquewords.txt'
xfull = np.loadtxt(wordsfp, dtype=int)
