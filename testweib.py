import numpy as np
import scipy.special as sp
import scipy.optimize as op
import fit
from scipy.stats import norm
import integration_constants as ic


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

"""
    Version 1: Binned
    f,F are for cts x, where x>=0
    g,G are for discrete x, where x>=0
    F(x; mu, sigma) = (1/2) (1+erf((log(x)-mu)/(sqrt(2)*sigma)))
    g(x; mu, sigma) = F(x) - F(x-1), x>=0
    G(x, mu, sigma) = F(x) = F(x) (because sum(g) is telescoping)
    h,H are for discrete, where x>= xmin
    h(x; mu, sigma) = (1/(0.5-F(xmin-1)))*g(x), x>=xmin
"""

def logpdf(x, a, b):
    xmin = np.min(x)
    F = lambda x: np.exp(-(x/b)**a)
    g = lambda x: F(x)-F(x+1)
    h = -np.log(F(xmin))+np.log(g(x))
    return h



terr = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/terrorism.txt'
words = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/uniquewords.txt'
xfull = np.loadtxt(words, dtype=int)
[alpha,xmin, ntail, Lpl, ks] = fit.pl(xfull)
x = xfull[xfull>=xmin]

# initial estimates
theta0 = initialguessweib(x)
# optimize
negloglike = lambda theta: -np.sum(logpdf(x,theta[0],theta[1]))
tol = 1E-5
bnds=[(tol,1),(0.01,None)]
res = op.minimize(negloglike, theta0, bounds=bnds, method='L-BFGS-B')
print(res)
theta = res.x
loglike = -negloglike(theta)
LplV = pllogpdf(x,alpha)
LstrexpV = logpdf(x,theta[0], theta[1])
R, p, normR = vuong(LplV, LstrexpV)
# check if statistically significant
decisionthresh=0.1
dexp = decide(R, p, decisionthresh)
print dexp, R, p, loglike, normR
