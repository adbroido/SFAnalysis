import sortgmls as sg
import numpy as np
import pandas as pd
import igraph
import os
import time

def listgmls(gmldirpath):
    fpV = []
    for root, dirs, files in os.walk(gmldirpath):
        # avoid gmls that are known to crash igraph.read()
        if 'bad' in dirs:
            dirs.remove('bad')
        if 'crash' in dirs:
            dirs.remove('crash')
        # avoid any other directories
        for avdir in avoiddirs:
            if avdir in dirs:
                dirs.remove(avdir)
        for name in files:
            # leave out the bipartite projections so we can make our own
            if name.endswith('.gml') and '1mode' not in name:
                fpV.append(os.path.join(root, name))
    return fpV

def gmlnames(gml_fpV):
    fnV = []
    for fp in gml_fpV:
        fnV.append(fp.split('/')[-1])
    return fnV


squal_fp = '/Volumes/Durodon/squalodon_gmls/'
clean_fp = '/Volumes/Durodon/clean_gmls/'
dur_fp = '/Volumes/Durodon/durodon_gmls/'
avoiddirs = []

squal_gmls = listgmls(squal_fp)
dur_gmls = listgmls(dur_fp)
clean_gmls = listgmls(clean_fp)

squal_gml_names = gmlnames(squal_gmls)
dur_gml_names = gmlnames(dur_gmls)
clean_gml_names = gmlnames(clean_gmls)

# find files that are in durodon but not in squal
set(dur_gml_names).difference(clean_gml_names)
