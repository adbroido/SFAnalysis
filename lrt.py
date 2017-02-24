import numpy as np
import scipy.optimize as op
import integration_constants as ic
import fit


def nonnested(x, Lpl):
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
    dexp = exp(x,Lpl)
    # compare log normal
    dln = ln(x,Lpl)
    # compare stretched exponential
    dstrexp = strexp(x,Lpl)

    return [dexp, dln, dstrexp]


def exp(x,Lpl):
    """
    Perform likelihood ratio test for exponetial distribution. First fits an
    exponential distribution to the data. A Vuong statistic is calculated from
    the likelihoods of the exponential fit and the power law fit. The
    convergence criteria for exponential fitting is met by either hitting the
    maximum number of iterations, or by reaching a likelihood value at which a
    decision is statistically possible.

    Input:
        x       ndarray, data set to be fit
        Lpl     float, log likelihood of power law fit

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
    res = op.minimize(negloglike,lam0)
    lam = np.asscalar(res.x)
    loglike = -negloglike(lam)
    # perform lrt: log-likelihood ratio between discrete power law
    # andexponential distribution
    
