import sortgmls as sg
import numpy as np
import pandas as pd
import igraph
import os
import time



def buildGMLcatalog(gmldirpath, df, avoiddirs, overwrite):
    """ Walks through the subdirectories of a root to find all gml files, then
    catalogs the relevant information about the contained networks.

    Input:
        gmldirpath              string, path to the root directory where gmls are
        df                      DataFrame, catalog of the existing gml files
        avoiddirs               list, contains names of any directories to avoid
                                    i.e. 'n7'
        overwrite               boolean, if true, forces overwrite of rows that
                                    already exist in the catalog. Otherwise only
                                    new files are added.



    Output:
        df                      DataFrame, catalog of the existing gml files

    """
    # make list of file paths to gmls
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
    # check which gmls are already in df
    if overwrite == False:
        newfpV = set(fpV).difference(set(df['fp_gml']))
    else:
        newfpV = fpV
    # update the catalog
    for fp in newfpV:
        g = igraph.read(fp)
        splitfp = fp.split('/')
        name = splitfp[-1]
        df.loc[name] = np.nan
        df.loc[name]['fp_gml'] = fp
        df.loc[name]['Domain'] = splitfp[-4]
        df.loc[name]['Subdomain'] = splitfp[-3]
        # cut out the 'n'
        df.loc[name]['Graph_order'] = int(splitfp[-2][1:])
        df.loc[name]['Weighted'] = sg.weighted(g, fp)
        df.loc[name]['Directed'] = sg.directed(g)
        df.loc[name]['Bipartite'] = sg.bipartite(g, fp)
        df.loc[name]['Multigraph'] = sg.multigraph(g)
        df.loc[name]['Multiplex'] = sg.multiplex(g)
        if (df.loc[name] == 'error').any():
            df = df.drop(name)
    return df



def preprocess(gmldirpath, init=False, avoiddirs=[], overwrite=False, timeit=False, save=True):
    """
    Sort the GMLs
    """
    if init:
        df = pd.DataFrame(columns=['Domain', 'Subdomain', 'Graph_order',
                'fp_gml', 'Weighted', 'Directed', 'Bipartite', 'Multigraph', 'Multiplex'])
    else:
        df = pd.read_pickle('gmlcatalog.p')
    start = time.time()
    df = buildGMLcatalog(gmldirpath, df, avoiddirs, overwrite)
    elapsed = time.time()-start
    if timeit:
        print 'Elapsed time is %s seconds' %elapsed
    if save:
        df.to_pickle('gmlcatalog.p')


if __name__ == '__main__':
    # run this:
    localpath = '/Users/annabroido/CU/gmls'
    externalpath = '/Volumes/Durodon/gmls'
    avoiddirs = ['n4', 'n5', 'n6', 'n7', 'n8', 'n9']
    avoiddirs=[]
    #df = buildGMLcatalog(externalpath, df, avoiddirs, overwrite=True)
    preprocess(externalpath, timeit=True, avoiddirs=avoiddirs)
