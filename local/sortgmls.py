import numpy as np
import igraph
import glob
import os
import pickle
import pandas as pd


# filepath to error file
errorfp = 'catalogerror.txt'

def weighted(g, fp=''):
    df_entry = 0
    if 'weight' in g.es.attributes():
        if len(np.unique(g.es['weight'])) >1:
            df_entry = 1
    elif 'value' in g.es.attributes():
        errormessage = "%s is weighted but has attribute 'value' instead of 'weight'\n" %fp
        known = False
        f = open(errorfp, 'a+')
        f.seek(0)
        for line in f:
            if line == errormessage:
                known = True
        if not known:
            f.write(errormessage)
        f.close()
        if len(np.unique(g.es['value'])) >1:
            df_entry = 1
    return df_entry

def multigraph(g):
    """ Check whether the graph g is a multigraph. At this step, multiplex
    graphs of any kind are also included. These will be separated later.
    """
    if g.has_multiple():
        df_entry = 1
    else:
        df_entry = 0
    return df_entry

def directed(g):
    if g.is_directed():
        df_entry = 1
    else:
        df_entry = 0
    # if it was already there, we just return it
    return df_entry

def multiplex(g):
    df_entry = 0
    # check if edges have an attribute other than weight or value
    attributes = set(g.es.attributes())
    weightattributes = set(['weight', 'value'])
    setdiff = attributes.difference(weightattributes)
    if len(setdiff)>0:
        df_entry = 1
    return df_entry

def bipartite(g, fp=None):
    if g.is_bipartite():
        if 'type' in g.vs.attributes():
            if len(set(g.vs['type'])) > 1:
                df_entry = 1
            else:
                df_entry = 0
        else:
            if not fp:
                errormessage = "%s is bipartite and has no attribute 'type'\n" %fp
                known = False
                f = open(errorfp, 'a+')
                f.seek(0)
                for line in f:
                    if line == errormessage:
                        known = True
                if not known:
                    f.write(errormessage)
                f.close()
                df_entry = 'error'
            else:
                df_entry = 0
    else:
        df_entry = 0
    return df_entry
