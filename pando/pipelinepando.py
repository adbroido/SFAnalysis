import sys
sys.path.insert(0, '/Users/anbr3575/LRTAnalysis/analysis/')
import numpy as np
import pandas as pd
import fit
import lrt
import importfiles as im
import sys

"""
Takes as input one degree sequence file.
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

def writelrterror(errormessage):
    errorfp = 'lrterror.txt'
    known = False
    f = open(errorfp, 'a+')
    f.seek(0)
    for line in f:
        if line == errormessage:
            known = True
    if not known:
        f.write(errormessage)
    f.close()

def runanalysis(fp):
    x = im.readdata(fp)
    # pull out the name for reference
    splitstr = fp.split('/')
    fn = splitstr[-1]
    analysis = pd.DataFrame(columns=['n', 'alpha', 'xmin', 'ntail', 'Lpl', 'ppl',
                                     'dexp', 'dln', 'dstrexp', 'dplwc'],
                                     index=[fn])
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
        if np.isnan(analysis.loc[fn]['alpha']):
            analysis.loc[fn]['n'] = n
            [alpha, xmin, ntail,  L, ks] = fit.pl(x)
            p = fit.plpval(x,alpha, xmin, ks)
            analysis.loc[fn]['alpha'] = alpha
            analysis.loc[fn]['xmin'] = xmin
            analysis.loc[fn]['ntail'] = ntail
            analysis.loc[fn]['Lpl'] = L
            analysis.loc[fn]['ppl'] = p
        if np.isnan(analysis.loc[fn]['dexp']):
            # compare the alternative distributions
            xmin = analysis.loc[fn]['xmin']
            alpha = analysis.loc[fn]['alpha']
            x = x[x>=xmin]
            # compare the non-nested alternatives, return the decisions for each
            decisionthresh = 0.1
            [dexp, dln, dstrexp] = lrt.nonnested(x,alpha, decisionthresh)
            if dexp == 2:
                errormessage = "Exponential didn't converge for %s \n" %fp
                writelrterror(errormessage)
            if dln == 2:
                errormessage = "Log-normal didn't converge for %s \n" %fp
                writelrterror(errormessage)
            if dstrexp == 2:
                errormessage = "Stretched exponential didn't converge for %s \n" %fp
                writelrterror(errormessage)
            # fit the nested alternatives
            dplwc = lrt.nested(x, alpha, decisionthresh)
            if dplwc == 2:
                errormessage = "PLWC didn't converge for %s \n" %fp
                writelrterror(errormessage)
            # update dataframe
            analysis.loc[fn]['dexp'] = dexp
            analysis.loc[fn]['dln'] = dln
            analysis.loc[fn]['dstrexp'] = dstrexp
            analysis.loc[fn]['dplwc'] = dplwc
    return analysis


if __name__ == '__main__':
    userinput = sys.argv
    fp = userinput[1]
    outputfp = userinput[2]
    analysis = runanalysis(fp)
    analysis.to_csv(outputfp)
