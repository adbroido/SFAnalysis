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
