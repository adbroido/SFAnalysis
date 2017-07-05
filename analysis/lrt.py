import numpy as np
import scipy.optimize as op
import scipy.special as sp
import integration_constants as ic
import fit
from scipy.stats import norm, chi2


""" Contains functions used in likelihood ratio tests, comparing the power-law
fit with alternative distributions. All can be called directly, but nested() and
nonnested() are the only we use.

"""


def pllogpdf(x,alpha):
    """ Point-wise log-pdf of the power-law distribution. For use in computing
    point-wise log-likelihood ratios.

    Input:
        x               ndarray, data to be fit
        alpha           float, power-law exponent, comes from best-fit

    Output:
        logpdf          float, point-wise log pdf

    """
    logpdf = np.log(ic.plconst(np.array([alpha]), np.min(x))) -alpha*np.log(x)
    return logpdf

def vuong(LplV, LaltV):
    """ Vuong test for log-likelihood ratio tests. Computes the ratio and an
    associated p-value.

    Input:
        LplV            ndarray, power-law pointwise log-likelihood
        LaltV           ndarray, alternative pointwise log-likelihood

    Output:
        R               float, likelihood ratio
        p2              float, 2-tail p-value
        normR           float, normalized ratio (so we can assume standard
                            normal for calculating p-val)

    """
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
    return R, p2, normR

def decide(logratio, p, decisionthresh):
    """ Takes p-value and log-likelihood ratio and makes a decision whether the
    test is conclusive. If so, it indicates whether the power-law or alternative
    is favored.

    Input:
        logratio            float, normalized log-likelihood ratio
        p                   float, associated p-value
        decisionthresh      float, threshold for rejecting.
                                Default in paper is decisionthresh = 0.1

    Output:
        d                   int, decision from the LRT
                                Decisions are:
                                    1   -  power law better
                                    0   -  inconclusive
                                   -1   -  alternative dist better

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
    """ Takes p-value and log-likelihood ratio and makes a decision whether the
    test is conclusive. If so, it indicates whether the power-law or alternative
    is favored. This is for nested decision, which cannot come out in favor of
    the power-law. The only nested alternative distribution we use is power-law
    with exponential cutoff.

    Input:
        logratio            float, normalized log-likelihood ratio
        p                   float, associated p-value
        decisionthresh      float, threshold for rejecting.
                                Default in paper is decisionthresh = 0.1

    Output:
        d                   int, decision from the LRT
                                Decisions are:
                                    1   -  power law better
                                    0   -  inconclusive
                                   -1   -  alternative dist better

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
    the likelihoods of the exponential fit and the power law fit.

    Input:
        x                   ndarray, data set to be fit (has been truncated to start
                                at xmin)
        LplV                ndarray, pointwise likelihood values for power-law fit
        decisionthresh      float, threshold for rejecting.
                                Default in paper is decisionthresh = 0.1

    Output:
        dexp                int, decision about exponential distribution
    """
    # perform lrt: Log-likelihood ratio between discrete power law and
    # exponential distribution. This is done pointwise so that we can use
    # Vuong's statistic to estimate the variance in the ratio
    [lam, LexpV, convstatus] = fit.exp(x)
    if convstatus == True:
        R, p, normR = vuong(LplV, LexpV)
        # check if statistically significant
        dexp = decide(normR, p, decisionthresh)
    else:
        dexp = 2
    return dexp

def ln(x,LplV, decisionthresh):
    """
    Perform likelihood ratio test for log normal distribution. First
    fits a log normal distribution to the data. A Vuong statistic is
    calculated from the likelihoods of the exponential fit and the power law
    fit.

    Input:
        x                   ndarray, data set to be fit
        LplV                ndarray, pointwise loglikelihood for power-law fit
        decisionthresh      float, threshold for rejecting.
                                Default in paper is decisionthresh = 0.1

    Output:
        dln                 int, decision about log-normal distribution
    """
    [theta,LlnV, convstatus] = fit.ln(x)
    if convstatus == True:
        R, p, normR = vuong(LplV, LlnV)
        # check if statistically significant
        dln = decide(normR, p, decisionthresh)
    else:
        dln = 2
    return dln

def strexp(x,LplV, decisionthresh):
    """
    Perform likelihood ratio test for stretched exponetial (Weibull)
    distribution. First fits a stretched exponential distribution to the data,
    then calculates a Vuong statistic from the likelihoods of the exponential
    fit and the power law fit.

    Input:
        x                   ndarray, data set to be fit
        LplV                ndarray, pointwise loglikelihood for power-law fit
        decisionthresh      float, threshold for rejecting.
                                Default in paper is decisionthresh = 0.1

    Output:
        dstrexp             int, decision about exponential distribution
    """
    [theta, LstrexpV, convstatus] = fit.strexp(x)
    if convstatus == True:
        R, p, normR = vuong(LplV, LstrexpV)
        # check if statistically significant
        dstrexp = decide(normR, p, decisionthresh)
    else:
        dstrexp = 2
    return dstrexp

def nested(x, alpha, decisionthresh=0.1):
    """
    Perform likelihood ratio tests for alternative distributions that are in
    the power law family.

    Input:
        x                   ndarray, data set to be fit. Assumed to only have values above
                                the xmin value from the power law fit
        alpha               float, best fit power-law parameter
        decisionthresh      float, threshold for rejecting null hypothesis
                                Default is 0.1

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
        x                   ndarray, data set to be fit. Assumed to only have
                                values above the best fit xmin value (from PL)
        alpha               float, best fit power-law parameter
        decisionthresh      float, threshold for rejecting null hypothesis
                                        Default is 0.1

    Output:
        dexp                int, decision about exponential distribution
        dln                 int, decision about log normal distribution
        dstrexp             int, decision about stretched exp distribution
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
