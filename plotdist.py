import numpy as np

import matplotlib.pyplot as plt

# weibull
def weib(x,a,b):
    xmin = np.min(x)
    y = (a/b) * np.exp((xmin/b)**a)*(x/b)**(a-1)*np.exp(- (x/b)**a)
    return y


terrfp = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/terrorism.txt'
wordsfp = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/pli-R-v0.0.3-2007-07-25/datasets/uniquewords.txt'
xfull = np.loadtxt(terrfp, dtype=int)

xmin = 12
x = xfull[xfull>=xmin]


theta = np.array([  0.92748832,  36.42255823])

plt.plot(x,weib(x,theta[0],theta[1]))
plt.hist(x)
plt.show()
