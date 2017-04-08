import numpy as np
import pandas as pd
import fit
import lrt
import importfiles as im

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

def runanalysis(degfiles, analysis, overwrite=False):
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
            if analysis.loc[fn]['dexp']=='' or overwrite == True:
                # compare the alternative distributions
                xmin = analysis.loc[fn]['xmin']
                alpha = analysis.loc[fn]['alpha']
                x = x[x>=xmin]
                # compare the non-nested alternatives, return the decisions for each
                decisionthresh = 0.1
                [dexp, dln, dstrexp] = lrt.nonnested(x,alpha, decisionthresh)
                # fit the nested alternatives
                dplwc = lrt.nested(x, alpha, decisionthresh)
                # update dataframe
                analysis.loc[fn]['dexp'] = dexp
                analysis.loc[fn]['dln'] = dln
                analysis.loc[fn]['dstrexp'] = dstrexp
                analysis.loc[fn]['dplwc'] = dplwc
    return analysis



#run analysis
analysis = pd.read_pickle('analysis.p')
analysis = analysis.query('ppl>0.1 & ntail>50 & Graph_order<4')
degdir = '/Users/annabroido/Dropbox/Research/LRTAnalysis/degreesequences/'
degfiles = degdir+analysis.index
analysis = runanalysis(degfiles, analysis)
analysis.to_pickle('analysis.p')
