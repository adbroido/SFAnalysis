import glob
import numpy as np
import importfiles as imp
import fit
from matplotlib import pyplot as plt
from scipy.stats import rv_discrete
import integration_constants as ic
from mpl_toolkits.mplot3d import Axes3D


class discexp(rv_discrete):
    "Discrete exponential distribution"
    def __init__(self,xmin=0):
        super(discexp,self).__init__(a=xmin)
    def _pmf(self, x, lam):
        xmin = self.a
        out =  (1-np.exp(-lam))*np.exp(lam*xmin)*np.exp(-lam*x)
        return out

class discpl(rv_discrete):
    "Discrete power law distribution"
    def __init__(self,xmin=1):
        super(discpl,self).__init__(a=xmin)
    def _pmf(self, x, alpha):
        xmin = self.a
        const = ic.plconst(np.array([alpha]), xmin)
        out = const*(x**(-alpha))
        return out

class discplwc(rv_discrete):
    "Discrete power law with cutoff distribution"
    def __init__(self,xmin=1):
        super(discplwc,self).__init__(a=xmin)
    def _pmf(self, x, alpha, lam):
        xmin = self.a
        const = ic.plwcconst(alpha, lam, xmin)
        out = const*(x**(-alpha))*(np.exp(-lam*x))
        return out


xmin0 = 1
alpha0 = 1.5
lam0 = 0.32
size = 40

exp = discexp(xmin0)
expdata = exp.rvs(lam=lam0,size=size)

pl = discpl(xmin0)
pldata = pl.rvs(alpha=alpha0, size=size)

plwc = discplwc(xmin0)
plwcdata = plwc.rvs(alpha=alpha0, lam=lam0, size=size)


# loglikelihood of exponential fit
def explogpdf(x,lam,xmin):
    result = np.log(1-np.exp(-lam))+lam*xmin - lam*x
    return result
def exploglike(x,lam,xmin):
    return np.sum(explogpdf(x,lam, xmin))

lams = np.arange(0.01,2,0.01)
expexpLs = np.array([exploglike(expdata,lam, xmin0) for lam in lams])
plexpLs = np.array([exploglike(pldata,lam, xmin0) for lam in lams])
plwcexpLs = np.array([exploglike(plwcdata,lam, xmin0) for lam in lams])


# loglikelihood of pl fit
def pllogpdf(x,alpha, xmin):
    C = ic.plconst(np.array([alpha]),xmin)
    result = -alpha*np.log(x)+np.log(C)
    return result
def plloglike(x,alpha,xmin):
    return np.sum(pllogpdf(x,alpha, xmin))
alphas = np.arange(1.01,2.51,0.01)
expplLs = np.array([plloglike(expdata,al, xmin0) for al in alphas])
plplLs = np.array([plloglike(pldata,al, xmin0) for al in alphas])
plwcplLs = np.array([plloglike(plwcdata,al, xmin0) for al in alphas])


# loglikelihood of plwc fit
def plwclogpdf(x,alpha,lam, xmin):
    C = ic.plwcconst(alpha,lam, xmin)
    result = np.log(C)-alpha*np.log(x)-lam*x
    return result
def plwcloglike(x,alpha,lam,xmin):
    return np.sum(plwclogpdf(x,alpha,lam, xmin))
alphas = np.arange(1.01,2.51,0.01)
lams = np.arange(0.01,1.01,0.01)
a,l = np.meshgrid(alphas,lams)
plwcparams = np.dstack((a,l)).swapaxes(0,1)
expplwcLs = np.array([[plwcloglike(expdata,al,lam, xmin0) for lam in lams] for al in alphas])
plplwcLs = np.array([[plwcloglike(pldata,al,lam, xmin0) for lam in lams] for al in alphas])
plwcplwcLs = np.array([[plwcloglike(plwcdata,al,lam, xmin0) for lam in lams] for al in alphas])


lamfit_exp = lams[expexpLs.argmax()]
lamfit_pl = lams[plexpLs.argmax()]
lamfit_plwc = lams[plwcexpLs.argmax()]
alfit_exp = alphas[expplLs.argmax()]
alfit_pl = alphas[plplLs.argmax()]
alfit_plwc = alphas[plwcplLs.argmax()]
# plwcfit_exp = alphas[expplwcLs.argmax()]
# plwcfit_pl = alphas[plplwcLs.argmax()]
# plwcfit_plwc = alphas[plwcplwcLs.argmax()]


print('estimated lambda for exponential data = %s' %lamfit_exp)
print('estimated lambda for pl data = %s' %lamfit_pl)
print('estimated lambda for plwc data = %s' %lamfit_plwc)


print('estimated alpha for exponential data = %s' %alfit_exp)
print('estimated alpha for pl data = %s' %alfit_pl)
print('estimated alpha for plwc data = %s' %alfit_plwc)

# print('estimated [alpha, lam] for exponential data = %s' %plwcfit_exp)
# print('estimated [alpha, lam] for pl data = %s' %plwcfit_pl)
# print('estimated [alpha, lam] for plwc data = %s' %plwcfit_plwc)



# plot
fig = plt.figure()
expplt = fig.add_subplot(1,2,1)
expplt.plot_surface(alphas,lams, expplwcLs)
expplt.set_title('PLWC fit to exp data, alpha = %s, lam = %s' %(alpha0, lam0))
plplt = fig.add_subplot(1,2,2)
plplt.plot(alphas,lams, plplwcLs)
plplt.set_title('PLWC fit to PL data')
plt.show()

# surf plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X,Y = np.meshgrid(alphas,lams)
Z = plplwcLs.reshape(Y.shape)
ax.plot_surface(A,L,Z)
plt.show()
