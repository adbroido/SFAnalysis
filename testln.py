import numpy as np
import scipy.special as sp
import scipy.optimize as op
import fit

"""
    Version 1: Binned
    f,F are for cts x, where x>0
    g,G are for discrete x, where x>=0
    F(x; mu, sigma) = (1/2) (1+erf((log(x)-mu)/(sqrt(2)*sigma)))
    g(x; mu, sigma) = F(x+1) - F(x), x>=0
    G(x, mu, sigma) = F(x+1)-F(0) = F(x+1) (because sum(g) is telescoping)
    h,H are for discrete, where x>= xmin
    h(x; mu, sigma) = (1/(1-F(xmin)))*g(x), x>=xmin
"""
def fitbin(x, meth):
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
    #theta0 = np.array([-110.93892217,   10.7925151 ])
    #theta0 = np.array([-48.2584935 ,   6.09574478]) #nelder-mead
    #theta0  = np.array([-32.49729189,   5.09234993]) # L-BFGS_B
    # optimize
    negloglike = lambda theta: -np.sum(logpdf(x,theta[0],theta[1]))
    tol = 1E-1
    bnds=[(-n/5,None),(tol,None)]
    if meth == 'Nelder-Mead':
        res = op.minimize(negloglike, theta0, method='Nelder-Mead')
    else:
        res = op.minimize(negloglike, theta0, bounds=bnds, method=meth)
    return [theta0, res]

def checkbin(xmin, theta):
    mu = theta[0]
    sigma = theta[1]
    xvec = np.arange(xmin ,xmin+100000)
    F = lambda x: (1+sp.erf((np.log(x)-mu)/(np.sqrt(2)*sigma)))/2
    g = lambda x: F(x)-F(x-1)
    h = (1/(1-F(xmin)))*g(xvec)
    tot = np.sum(h)
    return tot

"""
    Version 2: brute force summed
    f,F are for cts x, where x>=1
    g,G are for discrete x, where x>=xmin
    f(x) = (1/x) exp(-(log(x)-mu)^2/(sqrt(2)*sigma))
    C = (sum(f(x) from x = xmin to infinity )
    g(x) = (1/C)*f(x), x>=xmin

"""
# def fitsum(x, meth):
#     def logpdf(x, mu, sigma):
#         xmin = np.min(x)
#         f = lambda x: np.exp(-((np.log(x)-mu)**2)/(np.sqrt(2)*sigma))/x
#         xvec = np.arange(xmin, xmin+100000)
#         C = np.sum(f(xvec))
#         g = -np.log(C)+np.log(f(x))
#         return g
#     # initial estimates
#     mu0 = np.mean(np.log(x))
#     sigma0 = np.std(np.log(x))
#     theta0 = np.array([mu0, sigma0])
#     # optimize
#     negloglike = lambda theta: -np.sum(logpdf(x,theta[0],theta[1]))
#     tol = 1E-1
#     bnds=[(-20,None),(tol,None)]
#     if meth == 'Nelder-Mead':
#         res = op.minimize(negloglike, theta0, method='Nelder-Mead')
#     else:
#         res = op.minimize(negloglike, theta0, bounds=bnds, method=meth)
#     return [theta0, res]
#
# def checksum(xmin, theta):
#     mu = theta[0]
#     sigma = theta[1]
#     f = lambda x: np.exp(-((np.log(x)-mu)**2)/(np.sqrt(2)*sigma))/x
#     xvec = np.arange(xmin, xmin+100000)
#     C = np.sum(f(xvec))
#     g = (1/C)*f(xvec)
#     tot = np.sum(g)
#     return tot



terr = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/terrorism.txt'
words = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/uniquewords.txt'
xfull = np.loadtxt(words, dtype=int)
[alpha,xmin, ntail, Lpl, ks] = fit.pl(xfull)
x = xfull[xfull>=xmin]

# nelderbin = fitbin(x, 'Nelder-Mead')
# boundedbin = fitbin(x, 'L-BFGS-B')
# newtonbin = fitbin(x, 'TNC')
# print 'Nelder binned = %s' %nelderbin[1].x
# print 'Bounded binned = %s' %boundedbin[1].x
# print 'Newton binned = %s' %newtonbin[1].x

# neldersum = fitsum(x, 'Nelder-Mead')
# boundedsum = fitsum(x, 'L-BFGS-B')
# newtonsum = fitsum(x, 'TNC')
# print 'Nelder summed = %s' %neldersum[1].x
# print 'Bounded summed = %s' %boundedsum[1].x
# print 'Newton summed = %s' %newtonsum[1].x



# # PLOT TIME
# import matplotlib.pyplot as plt
#
# def pdfbin(x, theta):
#     mu = theta[0]
#     sigma = theta[1]
#     xmin = np.min(x)
#     F = lambda x: (1+sp.erf((np.log(x)-mu)/(np.sqrt(2)*sigma)))/2
#     g = lambda x: F(x)-F(x-1)
#     h = (1/(1-F(xmin-1)))*g(x)
#     return h
#
# def cdfbin(x, theta):
#     H = np.cumsum(pdfbin(x,theta))
#     return H
#
# def pdfsum(x, theta):
#     mu = theta[0]
#     sigma = theta[1]
#     xmin = np.min(x)
#     f = lambda x: np.exp(-((np.log(x)-mu)**2)/(np.sqrt(2)*sigma))/x
#     xvec = np.arange(xmin, xmin+100000)
#     C = np.sum(f(xvec))
#     g = (1/C)*f(x)
#     return g
#
# def cdfsum(x, theta):
#     H = np.cumsum(pdfsum(x,theta))
#     return H
#
# def pdfcts(x, theta):
#     mu = theta[0]
#     sigma = theta[1]
#     xmin = np.min(x)
#     f = lambda x: np.exp(-((np.log(x)-mu)**2)/(np.sqrt(2)*sigma))/x
#     C = np.sqrt(2/(np.pi*sigma**2))*(1/sp.erfc((np.log(xmin)-mu)/(np.sqrt(2)*sigma)))
#     g = C*f(x)
#     return g
#
# def cdfcts(x, theta):
#     H = np.cumsum(pdfcts(x,theta))
#     return H
#
# def plotpdf(x, theta):
#     # pdf plot
#     plt.plot(xvec, pdfbin(xvec,theta),label='Binned PMF')
#     plt.plot(xvec, pdfsum(xvec,theta), label='Summed PMF')
#     plt.plot(xvec, pdfcts(xvec,theta), label='Continuous PDF')
#     plt.legend()
#     plt.show()
#
# def plotcdf(x, theta):
#     # cdf plot
#     binnedcdf = plt.plot(xvec, cdfbin(xvec,theta), label='Binned CDF')
#     summedcdf = plt.plot(xvec, cdfsum(xvec,theta), label='Summed CDF')
#     #ctscdf = plt.plot(xvec, cdfcts(xvec,theta), label='Continuous CDF')
#     plt.legend()
#     plt.show()
#
# # MAKE PLOTS
# xvec = np.arange(xmin,xmin+10)
# theta = [7,3]
