import numpy as np
import scipy.optimize as op
import integration_constants as ic
import scipy.special as sp
import time

""" Contains functions used in fitting the power-law, exponential, log-normal,
Weibull (stretched exponential), and power-law with exponential cutoff, as well
as plpval() to find a p-value for the power-law fit. All should be called
directly. All distributions are discrete.

"""


def pl(x):
    """ Fits a tail-conditional power-law to a data set. This implements brute
    force optimization (grid search) instead of using a built in optimizer. The
    grid on alpha runs from alstart to alstart+shift. This is based on Aaron's
    plfit.m Matlab code (http://tuvalu.santafe.edu/~aaronc/powerlaws/).

    Input:
        x           ndarray, ndim = 1, dtype = integer


    Output:
        alpha        float, exponent on x, must be > 1
        xmin         int, starting point for power law tail, must be >= 1
        ntail        int, number of datapoints above (and including) xmin
        L            float, log likelihood of the returned fit
        ks           float, goodness of fit statistic (Kolmogorov-Smirnov)
    """
    # find the and sort unique possible xmin values
    xminV = np.trim_zeros(np.unique(x))
    # initialize array of the fits for every xmin
    fitV = np.zeros([len(xminV),2])

    start_time = time.time()
    # initialize vector of constants
    # where the xmins start
    xminprev = min(xminV) - 1
    # initialize array of possible alpha values
    alstart = 1.01
    shift = 9.50
    alphaV = np.arange(alstart,alstart+shift,0.01)
    zetaV = sp.zeta(alphaV)
    constV = zetaV
    # shift up to start at the smallest xmin
    for j in range(xminprev):
        constV += -(1+j)**(-alphaV)

    # loop over the xmin values at find the best fit at each
    for i in range(len(xminV)):
        xmin = xminV[i]
        xtail = x[x>=xmin]
        ntail = len(xtail)
        # optimize over alpha
        # find the corresponding array of conditional log likelihoods
        Ls = -alphaV*np.sum(np.log(xtail)) - ntail*np.log(constV)
        # pull out the location of the best fit alpha
        aind = Ls.argmax()
        # find what alpha value is at this index
        alpha = alphaV[aind]
        # compute the KS statistic
        # theoretical cdf
        cdf = np.cumsum(range(np.min(xtail), np.max(xtail)+1)**(-alpha)/constV[aind])
        #  binned data
        xhist = np.histogram(xtail,range(np.min(xtail), np.max(xtail)+2))
        # empirical cdf
        edf = np.cumsum(xhist[0])/float(ntail)
        # KS stat
        ks = np.max(np.abs(cdf-edf))
        # add this KS stat and alpha to the array of fits
        fitV[i] = np.array([ks, alpha])
        # update the constants
        for j in range(xmin-xminprev):
            constV += -(xmin+j)**(-alphaV)
        xminprev = xmin


    # pull out the index of the smallest KS stat
    ksind = fitV[:,0].argmin()
    ks = fitV[ksind,0]
    # find the corresponding xmin
    xmin = xminV[ksind]
    # find the corresponding alpha
    alpha = fitV[ksind,1]
    # evaluate the likelihood here
    xtail = x[x>=xmin]
    ntail = len(xtail)
    start_time = time.time()
    const = sp.zeta(alpha) - np.sum(np.arange(1,xmin)**(-alpha))
    L = -alpha * np.sum(np.log(xtail)) - ntail*np.log(const)
    # print "-------%s seconds -----------" %(time.time()-start_time)
    # print "alpha = %s" %alpha
    # print "xmin = %s" %xmin
    return [alpha,xmin, ntail, L, ks]

def plpval(x, alpha, xmin, gof):
    """ Finds p-value for the power-law fit using a KS test. This is based on
    Aaron's plpva.m Matlab code (http://tuvalu.santafe.edu/~aaronc/powerlaws/).

    Input:
        x            ndarray, ndim = 1, dtype = integer
        alpha        float, exponent on x, must be > 1
        xmin         int, starting point for power law tail, must be >= 1
        gof           float, goodness of fit statistic (Kolmogorov-Smirnov)


    Output:
        p            p-value of the returned fit (reject PL hypothesis for p<0.1)
    """
    # set desired precision level in p-value
    eps = 0.01
    #num_resamps = int(np.ceil((1./4)*eps**(-2)))
    num_resamps = 1000
    bootstraps = np.zeros(num_resamps)
    n = len(x)
    xmax = np.max(x)
    tailinds = x>=xmin
    xtail = x[tailinds]
    xhead = x[~tailinds]
    ntail = len(xtail)
    nhead = len(xhead)
    ptail = float(ntail)/n
    mmax = 20*xmax
    # set the tail of the pdf
    #const_tail = ic.plconst(np.array(alpha),xmin)
    const_tail = sp.zeta(alpha) - np.sum(np.arange(1,xmin)**(-alpha))
    pdf_tail = np.arange(xmin,mmax+1)**(-alpha)/const_tail # i.e.; end at mmax
    # pad this with zeros (rewrite if we don't need to do this)
    pdf = np.zeros(mmax+1)
    pdf[xmin:] = pdf_tail
    # clean up in case this is a huge array
    del pdf_tail
    # set the cdf. rows are x-val or cdf(xval). So cdf(x=10) is cdf[1,10]
    cdf = np.array( [ np.arange(mmax+1), np.cumsum(pdf) ] )
    # tack on a last entry
    cdf = np.concatenate( (cdf , np.array([[mmax+1,1]]).T) , axis = 1 )

    # semi-parametric bootstrap
    starttime = time.time()
    for resamp_ind in range(num_resamps):
        # non-parametric bootstrap from the head of x
        # count how many of n random numbers are in the head, based on the probability of being in the head of x
        nnewhead = n
        while nnewhead >= n:
            nnewhead = np.sum(np.random.rand(n)>ptail)
        headinds = np.array([np.floor(nhead*np.random.rand(nnewhead))],dtype=int)
        newhead = xhead[headinds][0]
        nnewtail  = n-nnewhead

        # parametric bootstrap for the powerlaw tail
        rtail = np.sort(np.random.rand(nnewtail))
        newtail = np.zeros(nnewtail, dtype = int)
        indrtail = 0
        indnewtail = 0
        for xval in range(xmin, mmax+2):
            while (indrtail < len(rtail)) and (rtail[indrtail] <= cdf[1, xval]):
                indrtail += 1
            newtail[indnewtail:indrtail] = xval
            indnewtail = indrtail
            if indnewtail > nnewtail:
                break
        # combine into new sample
        newx = np.concatenate((newhead, newtail))
        if (newx == np.zeros_like(newx)).all():
            import pdb; pdb.set_trace()
        # fit this new sample
        [newalpha, newxmin, newntail, newLpl, newgof] = pl(newx)
        # print where we are
        current_p = np.sum(bootstraps[0:resamp_ind]>=gof)/(float(resamp_ind+1))
        # print "[%s]    p = %f" %(resamp_ind, current_p)
        # store gof stat
        bootstraps[resamp_ind] = newgof
        # if it's taking forever and we can end, do it
        if time.time() - starttime > 500:
            if resamp_ind > num_resamps/20.:
                if current_p<0.05 or current_p>0.5:
                    print "current p = %s   elapsed time = %s" %(current_p, time.time()-starttime)
                    return current_p
    p = np.sum(bootstraps>=gof)/float(num_resamps)
    print "p = %.3f   elapsed time = %s" %(p, time.time()-starttime)
    return p

def exp(x):
    """ Fits a tail-conditional exponential to a data set. The data is assumed
    to begin at xmin. The logpdf is what is calculated and returned, as this is
    more relevant for likelihood calculations.

    Input:
        x            ndarray, ndim = 1, dtype = integer

    Output:
        lam          float, exponential rate, must be > 0
        LV           ndarray, pointwise log likelihood of the returned fit
        convstatus   Boolean, True if the fit converged, false if not
    """
    xmin = np.min(x)
    ntail = len(x)
    # define log pdf
    def logpdf(x,lam):
        result = np.log(1-np.exp(-lam))+lam*xmin - lam*x
        return result
    if len(np.unique(x))<2:
        # this means every value is equal to xmin
        # return dummy answers and say we don't converge
        lam = 0
        LV =0
        convstatus = False
    else:
        # Moment based estimate for optimzation
        lam0 = np.log(1+float(ntail)/np.sum(x-xmin))
        # define negative log likelihood, the function we wish to minimize
        negloglike = lambda lam: -np.sum(logpdf(x,lam))
        tol = 1E-9
        res = op.minimize(negloglike,lam0, bounds=[(tol,None)],method='L-BFGS-B')
        lam = np.asscalar(res.x)
        convstatus = res.success
        LV = logpdf(x,lam)
    return [lam, LV, convstatus]

def ln(x):
    """ Fits a tail-conditional log normal distribution to a data set.
    The data is assumed to begin at xmin. The logpdf is what is calculated and
    returned, as this is more relevant for likelihood calculations.
    Discretization is done by binning the continuous distrbution
    (see text for details)

    Input:
        x               ndarray, ndim = 1, dtype = integer

    Output:
        theta           ndarray, [mu, sigma] where mu is a float, the mean of the
                            distribution, unbounded (though we impose a bound), and sigma is a
                            float, the standard deviation, and must be > 0
        LV              ndarray, pointwise log likelihood of the returned fit
        convstatus      Boolean, True if the fit converged, false if not
    """
    xmin = np.min(x)
    ntail = len(x)
    # define log pdf
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
    # optimize
    negloglike = lambda theta: -np.sum(logpdf(x,theta[0],theta[1]))
    tol = 1E-1
    bnds=[(-n/5,None),(tol,None)]
    res = op.minimize(negloglike, theta0, bounds=bnds, method='L-BFGS-B')
    theta = res.x
    convstatus = res.success
    LV = logpdf(x,theta[0], theta[1])
    return [theta, LV, convstatus]

def plwc(x, alpha0=None):
    """ Fits a tail-conditional power-law with exponential cutoff to a data set.
    The data is assumed to begin at xmin. The logpdf is what is calculated and
    returned, as this is more relevant for likelihood calculations.

    Input:
        x           ndarray, ndim = 1, dtype = integer
        alpha0      float, power-law exponent (optional input)

    Output:
        alpha        float, exponent on x, must be > -1
        lam          float, exponential rate, must be > 0
        LV           ndarray, pointwise log likelihood of the returned fit
        convstatus   Boolean, True if the fit converged, false if not
    """
    xmin = np.min(x)
    ntail = len(x)
    # define log pdf
    def logpdf(x,alpha, lam):
        xmin = np.min(x)
        C = ic.plwcconst(alpha,lam, xmin)
        result = -np.log(C) - alpha*np.log(x) - lam*x
        return result
    # Estimates for optimzation
    if alpha0 is None:
        alpha0 = pl(x)[0]
    lam0 = exp(x)[0]
    theta0 = np.array([alpha0,lam0])
    # define negative log likelihood, the function we wish to minimize
    negloglike = lambda theta: -np.sum(logpdf(x,theta[0], theta[1]))
    tol = 1E-5
    bnds=[(-1+tol,None),(tol,None)]
    res = op.minimize(negloglike, theta0, bounds=bnds)
    # res = op.minimize(negloglike,theta0, method='Nelder-Mead')
    theta = res.x
    convstatus = res.success
    alpha = theta[0]
    lam = theta[1]
    LV = logpdf(x,alpha, lam)
    return [alpha, lam, LV, convstatus]

def strexp(x):
    """ Fits a tail-conditional stretched exponential distribution to a data set.
    The data is assumed to begin at xmin. The logpdf is what is calculated and
    returned, as this is more relevant for likelihood calculations.
    Discretization is done by binning the continuous distrbution
    (see text for details)

    Input:
        x           ndarray, ndim = 1, dtype = integer

    Output:
        theta           ndarray, [a,b], dtype=float. a (0<a<1) is the mean of the
                            distribution and b is the standard deviation (b > 0)
        LV              ndarray, pointwise log likelihood of the returned fit
        convstatus      Boolean, True if the fit converged, false if not
    """
    xmin = np.min(x)
    ntail = len(x)
    # define log pdf
    def initialguessweib(x):
        """
        The follow is not an efficient estimator of the shape, but serves to start
        the approximation process off.  (See Johnson and Katz, ch. 20, section 6.3,
        "Estimators Based on Distribution of log X".) This is taken from the R
        package by Cosma Shalizi at http://tuvalu.santafe.edu/~aaronc/powerlaws/).
        """
        xmin = np.min(x)
        n = len(x)
        shape = (np.sqrt(6)/np.pi)*np.std(np.log(x))
        scale = (np.sum(x**shape)/n)**(1/shape)
        return np.array([shape,scale])

    # define log pdf, for use in likelihood calculation
    def logpdf(x, a, b, tol_pad = 1E-8):
        xmin = np.min(x)
        F = lambda x: np.exp(-(x/b)**a)
        g = lambda x: F(x)-F(x+1)
        h = -np.log(F(xmin)+tol_pad)+np.log(g(x)+tol_pad)
        return h
    # initial estimates
    # initial estimates
    theta0 = initialguessweib(x)
    # optimize
    negloglike = lambda theta: -np.sum(logpdf(x,theta[0],theta[1]))
    tol = 1E-5
    bnds=[(tol,1),(0.01,None)]
    res = op.minimize(negloglike, theta0, bounds=bnds, method='L-BFGS-B')
    theta = res.x
    convstatus = res.success
    LV = logpdf(x,theta[0], theta[1])
    return [theta, LV, convstatus]
