import numpy as np
import igraph
import glob
import os
import pickle
import pandas as pd

""" Assorted functions to check whether a graph (as an igraph object) has
certain properties. All are meant to be called directly.

"""


# filepath to error file
errorfp = 'gmlerror.txt'

def weighted(g, fp=''):
    """ Check whether the graph g is weighted. The built-in igraph check only
    looks for a type label of 'weight'. Sometimes gmls will have 'value' instead
    so this makes sure to check for both. If 'value' is used, this gml is noted
    in an error file so we can fix it later.

    Input:
        g               igraph object, graph to be checked
        fp              string, filepath to gml file. To be used in error file
                        if necessary.

    Output:
        df_entry        int, 0 means not multiplex, 1 means multiplex

    """
    df_entry = 0
    if 'weight' in g.es.attributes():
        if len(np.unique(g.es['weight'])) >1:
            df_entry = 1
    elif 'value' in g.es.attributes():
        errormessage = "%s is weighted but has attribute 'value' instead of 'weight'\n" %fp
        # only add this line to error file if we  haven't already noted this
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

    Input:
        g               igraph object, graph to be checked

    Output:
        df_entry        int, 0 means not multigraph, 1 means multigraph

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
    """ Check whether the graph g is multiplex by checking for edge types.

    Input:
        g               igraph object, graph to be checked

    Output:
        df_entry        int, 0 means not multiplex, 1 means multiplex

    """
    df_entry = 0
    # check if edges have an attribute other than weight or value
    attributes = set(g.es.attributes())
    weightattributes = set(['weight', 'value'])
    setdiff = attributes.difference(weightattributes)
    if len(setdiff)>0:
        df_entry = 1
    return df_entry

def bipartite(g, fp=None):
    """ Check whether the graph g is bipartite.

    Input:
        g               igraph object, graph to be checked
        fp              string, path to gml file

    Output:
        df_entry        int or string, 0 means not bipartite, 1 means bipartite
                        'error' means the gml file is not structured correctly

    """
    if g.is_bipartite():
        if 'type' in g.vs.attributes():
            if len(set(g.vs['type'])) > 1:
                df_entry = 1
            else:
                df_entry = 0
        else:
            if fp:
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
