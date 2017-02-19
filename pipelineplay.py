import numpy as np
import pandas as pd
import fit
import lrt
import importfiles as im
import igraph
import sortgmls as sg

"""
Assumes  we've already got the degree sequence files.
"""

def writeerror(errormessage):
    errorfp = 'analysiserror.txt'
    known = False
    f = open(errorfp, 'a+')
    f.seek(0)
    for line in f:
        if line == errormessage:
            known = True
    if not known:
        f.write(errormessage)
    f.close()

def powerlawanalysis(degfiles, analysis, overwrite=False):
    for fp in degfiles:
        x = im.readdata(fp)
        # pull out the name for reference
        splitstr = fp.split('/')
        fn = splitstr[-1]
        # catch if it's not in there
        if fn not in analysis.index:
            print 'this degree sequence has not been seen before'
        # note if there is a problem with the file
        n = len(x)
        if np.mean(x) < 2 or np.mean(x) > np.sqrt(n):
            errormessage = "%s has a bad mean degree \n" %fp
            writeerror(errormessage)
        # catch for trivial degree sequences
        elif len(np.unique(x)) == 1:
            errormessage = "%s contains only one unique value \n" %fp
            writeerror(errormessage)
        else:
            if analysis.loc[fn]['alpha']=='' or overwrite == True:
                analysis.loc[fn]['n'] = n
                [alpha, xmin, ntail,  L, ks] = fit.pl(x)
                p = fit.plpval(x,alpha, xmin, ks)
                analysis.loc[fn]['alpha'] = alpha
                analysis.loc[fn]['xmin'] = xmin
                analysis.loc[fn]['ntail'] = ntail
                analysis.loc[fn]['LPL'] = L
                analysis.loc[fn]['pPL'] = p
    return analysis


def lrtanalysis(degfiles, analysis, overwrite=False):
    for fp in degfiles:
        x = im.readdata(fp)
        # pull out the name for reference
        splitstr = fp.split('/')
        fn = splitstr[-1]
        # catch if it's not in there
        if fn not in analysis.index:
            print 'this degree sequence has not been seen before'
        # note if there is a problem with the file
        if np.mean(x) < 2 or np.mean(x) > np.sqrt(n):
            errormessage = "%s has a bad mean degree \n" %fp
            # don't write again
        # catch for trivial degree sequences
        elif len(np.unique(x)) == 1:
            errormessage = "%s contains only one unique value \n" %fp
            # don't write again
        else:
            if analysis.loc[fn]['dexp']=='' or overwrite == True:
                xmin = analysis.loc[fn]['xmin']
                # truncate the data
                x = x[x>=xmin]
                # compare the alternative distributions
                # compare the non-nested alternatives, return the decisions for each
                [dexp, dln, dstrexp] = lrt.nonnested(x, Lpl)
                # fit the nested alternatives
                [dplwc] = lrt.nested(x, Lpl)
                # update dataframe
                analysis.loc[fn]['dexp'] = dexp
                analysis.loc[fn]['dln'] = dln
                analysis.loc[fn]['dstrexp'] = dstrexp
                analysis.loc[fn]['dplwc'] = dplwc
    return analysis




#run PL analysis
analysis = pd.read_pickle('analysis.p')
degdir = '/Users/annabroido/Dropbox/Research/LRTAnalysis/degreesequences/'
degfiles = degdir+analysis.index
analysis = powerlawanalysis(degfiles, analysis)
analysis.to_pickle('analysis.p')

# now run the LRT analysis. Later make these one step
analysis = lrtanalysis(degfiles, analysis)
analysis.to_pickle('analysis.p')



# compare the alternative distributions
# compare the non-nested alternatives, return the decisions for each
[dexp, dln, dstrexp] = lrt.nonnested(x, Lpl)
# fit the nested alternatives
[dpls, dplwc] = lrt.nested(x, Lpl)
# save the output
analysis = pd.Series([name, xmin, dexp, dln, dstrexp, dpls, dplwc])
