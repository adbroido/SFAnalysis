import sys
sys.path.insert(0, '/Users/anbr3575/LRTAnalysis/analysis/')
import numpy as np
import pandas as pd
import fit
import lrt
import importfiles as im

"""
Takes as input one degree sequence file. Runs power-law fit, p-test, and LRTs,
and records results in a csv file.
"""
def writeerror(errormessage):
    """ Writes explanation of error if the mean degree is bad.
    """
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
    """ Writes explanation of error if any lrt test fails because the fit
    didnt converge for the alternative distribution.

    """
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

def runanalysis(fp, outputfp):
    """ Reads a degree sequence file and runs power-law and LRT analysis on it.
    Results are stored in a DataFrame, which is saved to a csv.

    Input:
        fp                      path to degree sequence file
        outputfp                path to results output file


    Output:
        analysis                output is written to a csv file, not returned

    """
    x = im.readdata(fp)
    # pull out the name for reference
    splitstr = fp.split('/')
    fn = splitstr[-1]
    analysis = pd.DataFrame(columns=['n', 'alpha', 'xmin', 'ntail', 'Lpl', 'ppl',
                                     'dexp', 'dln', 'dstrexp', 'dplwc'],
                                     index=[fn])
    # note if there is a problem with the degree sequence. Should have already been caught
    n = len(x)
    if np.mean(x) < 2 or np.mean(x) > np.sqrt(n):
        errormessage = "%s has a bad mean degree \n" %fp
        writeerror(errormessage)
    # catch for trivial degree sequences
    elif len(np.unique(x)) == 1:
        errormessage = "%s contains only one unique value \n" %fp
        writeerror(errormessage)
    else:
        analysis.loc[fn]['n'] = n
        # pl analysis
        [alpha, xmin, ntail,  L, ks] = fit.pl(x)
        p = fit.plpval(x,alpha, xmin, ks)
        analysis.loc[fn]['alpha'] = alpha
        analysis.loc[fn]['xmin'] = xmin
        analysis.loc[fn]['ntail'] = ntail
        analysis.loc[fn]['Lpl'] = L
        analysis.loc[fn]['ppl'] = p
        # compare the alternative distributions
        xmin = analysis.loc[fn]['xmin']
        alpha = analysis.loc[fn]['alpha']
        x = x[x>=xmin]
        # compare the non-nested alternatives, return the decisions for each
        decisionthresh = 0.1 # significance level to reject null hypothesis
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
    analysis.to_csv(outputfp)


if __name__ == '__main__':
    userinput = sys.argv
    # file path to degree sequence
    fp = userinput[1]
    # path to output csv file
    outputfp = userinput[2]
    runanalysis(fp, outputfp)
