import numpy as np
import igraph
import glob
import os
import pickle
import fit
import sortgmls as sg
import pandas as pd
import importfiles as im


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
        meandeg = np.mean(x)
        if meandeg < 2:
            errormessage = "%s mean degree is too small: %f \n" %(fp, meandeg)
            writeerror(errormessage)
        elif np.mean(x) > np.sqrt(n):
            errormessage = "%s's mean degree is too large: %f \n" %(fp, meandeg)
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
                analysis.loc[fn]['Lpl'] = L
                analysis.loc[fn]['ppl'] = p
    return analysis



analysis = pd.read_pickle('analysis.p')
degdir = '/Users/annabroido/Dropbox/Research/LRTAnalysis/degreesequences/'
degfiles = degdir+analysis.index
analysis = powerlawanalysis(degfiles, analysis)
analysis.to_pickle('analysis.p')
