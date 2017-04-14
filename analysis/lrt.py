import numpy as np
import scipy.optimize as op
import scipy.special as sp
import integration_constants as ic
import fit
from scipy.stats import norm, chi2


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
    """
    Decisions are   1  -   power law better
                    0  -   inconclusive
                   -1  -   alternative dist better
    """
    if p <= decisionthresh:
        if logratio == 0:
            d = 0
        elif logratio > 0:
            # PL better
            d = 1
        else:
            # alt better
            d = -1
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

def exp(x, LplV, decisionthresh):
    """
    Perform likelihood ratio test for exponetial distribution. First fits an
    exponential distribution to the data. A Vuong statistic is calculated from
    the likelihoods of the exponential fit and the power law fit. The
    convergence criteria for exponential fitting is met by either hitting the
    maximum number of iterations, or by reaching a likelihood value at which a
    decision is statistically possible.

    Input:
        x           ndarray, data set to be fit (has been truncated to start
                    at xmin)
        LplV        ndarray, pointwise likelihood values for power-law fit

    Output:
        dexp    int, decision about exponential distribution
    """
    # perform lrt: Log-likelihood ratio between discrete power law and
    # exponential distribution. This is done pointwise so that we can use
    # Vuong's statistic to estimate the variance in the ratio
    [lam, LexpV, convstatus] = fit.exp(x)
    if convstatus == True:
        R, p, normR = vuong(LplV, LexpV)
        # check if statistically significant
        dexp = decide(R, p, decisionthresh)
    else:
        dexp = 2
    return dexp

def ln(x,LplV, decisionthresh):
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
        dln    int, decision about log-normal distribution
    """
    [theta,LlnV, convstatus] = fit.ln(x)
    if convstatus == True:
        R, p, normR = vuong(LplV, LlnV)
        # check if statistically significant
        dln = decide(R, p, decisionthresh)
    else:
        dln = 2
    return dln

def strexp(x,LplV, decisionthresh):
    """
    Perform likelihood ratio test for stretched exponetial distribution. First
    fits a stretched exponential (weibull) distribution to the data. A Vuong
    statistic is calculated from the likelihoods of the exponential fit and the
    power law fit. The convergence criteria for exponential fitting is met by
    either hitting the maximum number of iterations, or by reaching a likelihood
    value at which a decision is statistically possible.

    Note: There are several functions which all claim to be the "discrete
    Weibull" distribution, this code uses the Nakagawa-Osaki discretization,
    see http://ljk.imag.fr/membres/Olivier.Gaudoin/ICRSA03Gaudoin.pdf for
    lists of others

    Input:
        x           ndarray, data set to be fit
        alpha       float, exponent of power law fit

    Output:
        dstrexp    int, decision about exponential distribution
    """
    [theta, LstrexpV, convstatus] = fit.strexp(x)
    if convstatus == True:
        R, p, normR = vuong(LplV, LstrexpV)
        # check if statistically significant
        dstrexp = decide(R, p, decisionthresh)
    else:
        dstrexp = 2
    return dstrexp

def nested(x, alpha, decisionthresh=0.1):
    """
    Perform likelihood ratio tests for alternative distributions that are not
    in the power law family.

    Input:
        x       ndarray, data set to be fit. Assumed to only have values above
                the xmin value from the power law fit

    Output:
        dplwc    int, decision about power law with exponential cutoff
                    Decisions are   1  -   power law better
                                    0  -   inconclusive
                                   -1  -   alternative dist better
    """
    LplV = pllogpdf(x,alpha)
    Lpl = np.sum(LplV)
    # compare plwc
    [alpha, lam, LplwcV, convstatus] = fit.plwc(x, alpha)
    if convstatus == True:
        Lplwc = np.sum(LplwcV)
        R = Lpl-Lplwc
        p = 1-chi2.cdf(-2*R, df=1)
        dplwc = decidenested(R, p, 0.1)
    else:
        dplwc = 2
    return dplwc

def nonnested(x, alpha, decisionthresh=0.1):
    """
    Perform likelihood ratio tests for alternative distributions that are not
    in the power law family.

    Input:
        x       ndarray, data set to be fit. Assumed to only have values above
                the xmin value from the power law fit

    Output:
        dexp    int, decision about exponential distribution
        dln     int, decision about log normal distribution
        dstrexp int, decision about stretched exponential distribution
                    Decisions are   1  -   power law better
                                    0  -   inconclusive
                                   -1  -   alternative dist better
    """
    LplV = pllogpdf(x,alpha)
    # compare exponential
    dexp = exp(x,LplV, decisionthresh)
    # compare log normal
    dln = ln(x,LplV, decisionthresh)
    # compare stretched exponential
    dstrexp = strexp(x,LplV, decisionthresh)

    return [dexp, dln, dstrexp]
