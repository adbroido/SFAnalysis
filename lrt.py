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
    res = op.minimize(negloglike,lam0, bounds=[(tol,None)],method='L-BFGS-B')
    lam = np.asscalar(res.x)
    loglike = -negloglike(lam)
    # perform lrt: Log-likelihood ratio between discrete power law and
    # exponential distribution. This is done pointwise so that we can use
    # Vuong's statistic to estimate the variance in the ratio
    LplV = pllogpdf(x,alpha)
    LexpV = logpdf(x,lam)
    R, p, normR = vuong(LplV, LexpV)
    # check if statistically significant
    dexp = decide(R, p, decisionthresh)
    return dexp, R, p, loglike, normR

def strexp(x,alpha, decisionthresh):
    """
    Perform likelihood ratio test for stretched exponetial distribution. First
    fits a stretched exponential distribution to the data. A Vuong statistic is
    calculated from the likelihoods of the exponential fit and the power law
    fit. The convergence criteria for exponential fitting is met by either
    hitting the maximum number of iterations, or by reaching a likelihood value
    at which a decision is statistically possible.

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
    # x has already been truncated to start at xmin
    # define cts pdf for discretization
    def initialguessweib(x):
        """
        The follow is not an efficient estimator of the shape, but serves to start
        the approximation process off.  (See Johnson and Katz, ch. 20, section 6.3,
        "Estimators Based on Distribution of log X".)
        """
        xmin = np.min(x)
        n = len(x)
        shape = (np.sqrt(6)/np.pi)*np.std(np.log(x))
        scale = (np.sum(x**shape)/n)**(1/shape)
        return np.array([shape,scale])

    def logctscdf(x,a,b):
        """
        Calculate the log cumulative distribution function (pdf) of data from a
        tail-conditional continuous stretched exponential:

            f(x) = (a/b) exp((xmin/b)^a) (x/b)^(a-1) exp(- (x/b)^a)
            F(x) = 1- exp((xmin/b)^a)*exp(-(x/b)^a)

        where a is shape and b is scale. This function returns the log of the
        upper tail cdf, ie P(X > x) = 1 - P(X>=x):
            F(x) = exp((xmin/b)^a)*exp(-(x/b)^a)
            lF(x)= (xmin/b)^a - (x/b)^a

        Inputs:
        x = vector of data to be fit. Assumed to be truncated (only x>=xmin)
        a = shape parameter
        b = scale paramter

        Output:
        lF = cdf of x
        """
        xmin = np.min(x)
        lF = (xmin/b)**a - (x/b)**a
        return lF

    # define log pdf, for use in likelihood calculation
    def logpdf(x, a, b):
        """
        compute PMF as increments of the continuous distribution function. All
        distributions are tail-conditional.
        Do calculations on a logarithmic scale:
        Let log(b) - log(a) = h = log(b/a), b > a
        Then b-a = a(b/a - 1) = a(exp(log(b/a)) - 1) = a(exp(h) - 1)
        log(b-a) = log(a) + log(exp(h) - 1)
        """
        lF1 = logctscdf(x,a,b)
        lF2 = logctscdf(x+1,a,b)
        result = lF2 + np.log(np.exp(lF1-lF2)-1)
        return result


    # initial parameter estimate for optimzation
    theta0 = initialguessweib(x)
    # define negative log likelihood, the function we wish to minimize
    negloglike = lambda theta: -np.sum(logpdf(x,theta[0],theta[1]))
    tol = 1E-5
    bounds=[(tol,None),(tol,None)]
    res = op.minimize(negloglike,theta0, method='Nelder-Mead')
    theta = res.x
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
    # x has already been truncated to start at xmin
    def arg(x):
        arg = (np.log(x)-mu)/(np.sqrt(2)*s)
        return arg
    # define log likelihood
    def lnloglike(x, mu, s):
        """
        compute PMF as increments of the continuous distribution function. All
        distributions are tail-conditional.
        """
        xmin = np.min(x)
        n = len(x)
        # 1/constant of normalization is P(X>=xmin), since discretized use P(X>=xmin-0.5)
        Cm = (1./2)*(1-sp.erf(arg(xmin-0.5)))
        L = -n*np.log(Cm) + np.sum(np.log(0.5*(sp.erfc(arg(x-0.5))-sp.erfc(arg(x+0.5)))))
        return L

    def logpdf(x,mu,s):
        Cm = (1./2)*(1-sp.erf(arg(xmin-0.5)))
        f = -np.log(Cm) + np.log(0.5*(sp.erfc(arg(x-0.5))-sp.erfc(arg(x+0.5))))
        return f

    # initial parameter estimate for optimzation
    theta0 = [np.mean(np.log(x)),np.std(np.log(x))]
    # define negative log likelihood, the function we wish to minimize
    negloglike =lambda theta: -loglikelihood(x,theta[0],theta[1])
    tol = 1E-5
    bounds=[(tol,None),(tol,None)]
    res = op.minimize(negloglike,theta0, method='Nelder-Mead')
    theta = np.asscalar(res.x)
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
