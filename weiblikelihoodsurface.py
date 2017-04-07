import numpy as np
import scipy.special as sp
import scipy.optimize as op
import fit

def loglikebin(x,theta):
    a = theta[0]
    b = theta[1]
    n = len(x)
    xmin = np.min(x)
    F = lambda x: np.exp(-(x/float(b))**a) # upper tail probability, ie 1-cdf
    g = lambda x: F(x)- F(x+1)
    L = -n*np.log(F(xmin)) + np.sum(np.log(g(x)))
    return L


def loglikesum(x,theta):
    a = theta[0]
    b = theta[1]
    n = len(x)
    xmin = np.min(x)
    f = lambda x: (float(a)/b)*(x/float(b))**(a-1)*np.exp(-(x/float(b))**a)
    xvec = np.arange(xmin, xmin+100000)
    C = np.sum(f(xvec))
    L = -n*np.log(C) + np.sum(np.log(f(x)))
    return L


def findmax(L):
Lmax = np.max(L)
inds = np.where(L==Lmax)
aind = inds[0][0]
bind = inds[1][0]
    return Lmax,aind,bind


terr = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/terrorism.txt'
words = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/uniquewords.txt'
xfull = np.loadtxt(words, dtype=int)
[alpha,xmin, ntail, Lpl, ks] = fit.pl(xfull)
x = xfull[xfull>=xmin]

# binned from 0
tol = 1E-1
amin = tol
amax = 1-tol
aV = np.linspace(amin, amax, 100)
bmin = 10
bmax = 2*xmin
bV = np.linspace(bmin, bmax, 100)
Lbin = np.zeros((len(aV), len(bV)))
for i in range(len(aV)):
    for j in range(len(bV)):
        theta = np.array([aV[i], bV[j]])
        Lbin[i,j] = loglikebin(x,theta)

# np.savetxt('weiblikebinwords.csv', Lbin, delimiter=',')

# summed from 0
Lsum = np.zeros((len(aV), len(bV)))
for i in range(len(aV)):
    for j in range(len(bV)):
        theta = np.array([aV[i], bV[j]])
        Lsum[i,j] = loglikesum(x,theta)
# np.savetxt('weiblikesum.csv', Lsum, delimiter=',')




# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# X = np.outer(muV, np.ones(np.size(sigmaV)))
# Y = np.outer(np.ones(np.size(muV)), sigmaV)
# ax.plot_surface(X, Y, L)
