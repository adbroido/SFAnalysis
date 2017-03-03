import numpy as np
import scipy.optimize as op
import integration_constants as ic
import fit
from scipy.stats import norm


def nonnested(x, Lpl, decisionthresh=0.1):
    """
    Perform likelihood ratio tests for alternative distributions that are not
    in the power law family.

    Input:
        x       ndarray, data set to be fit. Assumed to only have values above
                the xmin value from the power law fit
        Lpl     float, log likelihood of power law fit

    Output:
        dexp    int, decision about exponential distribution
        dln     int, decision about log normal distribution
        dstrexp int, decision about stretched exponential distribution
                    Decisions are  -1  -   power law better
                                    0  -   inconclusive
                                    1  -   alternative dist better
    """
    # compare exponential
    dexp = exp(x,Lpl, decisionthresh)
    # compare log normal
    dln = ln(x,Lpl, deicisionthresh)
    # compare stretched exponential
    dstrexp = strexp(x,Lpl, decisionthresh)

    return [dexp, dln, dstrexp]

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

def decide(logratio, p, decisionthresh):
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


def exp(x,alpha, decisionthresh):
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
    # x has already been truncated to start at xmin
    xmin = np.min(x)
    ntail = len(x)
    # define log pdf, for use in likelihood calculation
    def logpdf(x,lam):
        result = np.log(1-np.exp(-lam))+lam*xmin - lam*x
        return result
    # Moment based estimate for optimzation
    lam0 = np.log(1+float(ntail)/np.sum(x-xmin))
    # define negative log likelihood, the function we wish to minimize
    negloglike = lambda lam: -np.sum(logpdf(x,lam))
    tol = 1E-9
    res = op.minimize(negloglike,lam0, bounds=[(tol,None)])
    lam = np.asscalar(res.x)
    loglike = -negloglike(lam)
    # perform lrt: Log-likelihood ratio between discrete power law and
    # exponential distribution. This is done pointwise so that we can use
    # Vuong's statistic to estimate the variance in the ratio
    # NOTE: MIGHT NOT NEED SOME OF THIS ABOVE STUFF
    LplV = pllogpdf(x,alpha)
    LexpV = logpdf(x,lam)
    R, p, normR = vuong(LplV, LexpV)
    # check if statistically significant
    dexp = decide(R, p, decisionthresh)
    return dexp, R, p, loglike, normR

def exp(x,alpha, decisionthresh):
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
    # x has already been truncated to start at xmin
    xmin = np.min(x)
    ntail = len(x)
    # define log pdf, for use in likelihood calculation
    def logpdf(x,lam):
        result = np.log(1-np.exp(-lam))+lam*xmin - lam*x
        return result
    # Moment based estimate for optimzation
    lam0 = np.log(1+float(ntail)/np.sum(x-xmin))
    # define negative log likelihood, the function we wish to minimize
    negloglike = lambda lam: -np.sum(logpdf(x,lam))
    tol = 1E-9
    res = op.minimize(negloglike,lam0, bounds=[(tol,None)])
    lam = np.asscalar(res.x)
    loglike = -negloglike(lam)
    # perform lrt: Log-likelihood ratio between discrete power law and
    # exponential distribution. This is done pointwise so that we can use
    # Vuong's statistic to estimate the variance in the ratio
    # NOTE: MIGHT NOT NEED SOME OF THIS ABOVE STUFF
    LplV = pllogpdf(x,alpha)
    LexpV = logpdf(x,lam)
    R, p, normR = vuong(LplV, LexpV)
    # check if statistically significant
    dexp = decide(R, p, decisionthresh)
    return dexp, R, p, loglike, normR
