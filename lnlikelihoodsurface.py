import numpy as np
import scipy.special as sp
import scipy.optimize as op
import fit

def loglike0(x,theta):
    mu = theta[0]
    sigma = theta[1]
    n = len(x)
    xmin = np.min(x)
    F = lambda x: (sp.erfc((np.log(x)-mu)/(np.sqrt(2)*sigma)))/2
    g = lambda x: F(x)- F(x+1)
    L = -n*np.log(F(xmin)) + np.sum(np.log(g(x)))
    return L

def loglike1(x,theta):
    mu = theta[0]
    sigma = theta[1]
    n = len(x)
    xmin = np.min(x)
    F = lambda x: (sp.erfc((np.log(x)-mu)/(np.sqrt(2)*sigma)))/2
    g = lambda x: F(x-1)- F(x)
    L = -n*np.log(F(xmin-1)) + np.sum(np.log(g(x)))
    return L

def loglikesum(x,theta):
    mu = theta[0]
    sigma = theta[1]
    n = len(x)
    xmin = np.min(x)
    f = lambda x: np.exp(-((np.log(x)-mu)**2)/(np.sqrt(2)*sigma))/x
    xvec = np.arange(xmin, xmin+100000)
    C = np.sum(f(xvec))
    L = -n*np.log(C) + np.sum(np.log(f(x)))
    return L


def getloglike0(x):
    tol = 1E-1
    mumin = tol
    mumax = 10
    muV = np.linspace(mumin, mumax, 200)
    sigmamin = 1
    sigmamax = 10
    sigmaV = np.linspace(sigmamin, sigmamax, 100)
    L = np.zeros((len(muV), len(sigmaV)))
    for i in range(len(muV)):
        for j in range(len(sigmaV)):
            theta = np.array([muV[i], sigmaV[j]])
            L[i,j] = loglike0(x,theta)
    return L

def findmax(L):
    Lmax = np.max(L)
    inds = np.where(L==Lmax)
    muind = inds[0][0]
    sigind = inds[1][0]
    return Lmax,muind,sigind


terr = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/terrorism.txt'
words = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/uniquewords.txt'
xfull = np.loadtxt(terr, dtype=int)
[alpha,xmin, ntail, Lpl, ks] = fit.pl(xfull)
x = xfull[xfull>=xmin]

# binned from 0
tol = 1E-1
mumin = -20
mumax = 10
muV = np.linspace(mumin, mumax, 200)
sigmamin = 1
sigmamax = 10
sigmaV = np.linspace(sigmamin, sigmamax, 100)
L0 = np.zeros((len(muV), len(sigmaV)))
for i in range(len(muV)):
    for j in range(len(sigmaV)):
        theta = np.array([muV[i], sigmaV[j]])
        L0[i,j] = loglike0(x,theta)
# np.savetxt('lnlike0.csv', L0, delimiter=',')

# summed from 0
Lsum = np.zeros((len(muV), len(sigmaV)))
for i in range(len(muV)):
    for j in range(len(sigmaV)):
        theta = np.array([muV[i], sigmaV[j]])
        Lsum[i,j] = loglikesum(x,theta)
# np.savetxt('lnlikesum.csv', Lsum, delimiter=',')

# binned from 1
tol = 1E-1
mumin = tol
mumax = 10
muV = np.linspace(mumin, mumax, 200)
sigmamin = 1
sigmamax = 10
sigmaV = np.linspace(sigmamin, sigmamax, 100)
L1 = np.zeros((len(muV), len(sigmaV)))
for i in range(len(muV)):
    for j in range(len(sigmaV)):
        theta = np.array([muV[i], sigmaV[j]])
        L1[i,j] = loglike1(x,theta)
# np.savetxt('lnlike1.csv', L1, delimiter=',')
#
#
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# X = np.outer(muV, np.ones(np.size(sigmaV)))
# Y = np.outer(np.ones(np.size(muV)), sigmaV)
# ax.plot_surface(X, Y, L)
