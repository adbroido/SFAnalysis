import sortgmls as sg
import numpy as np
import pandas as pd
import igraph
import os
import time

""" Builds a catalog of gmls organized by type taking as input just the path to
all of the gml directories. preprocess() is the only function designed to be
called directly.

"""


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
        print fp
        g = igraph.read(fp)
        splitfp = fp.split('/')
        name = splitfp[-1]
        # add new row or overwrite existing row
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
            # this catches bad bipartite gmls
            df = df.drop(name)
    return df



def preprocess(gmldirpath, init=False, avoiddirs=[], overwrite=False, timeit=False, save=True):
    """
    Categorize the properties of the GMLs

    Input:
        gmldirpath                  string, path to outer gml directory
                                    (containing all Domain folders and subfolders)
        init                        Boolean, if True the catalog is rebuilt from
                                    scratch
        avoiddirs                   list, directories to avoid. For example 'n7'
                                    Useful for avoiding bigger directories for
                                    smaller preliminary analysis
        overwrite                   Boolean, if True, all rows will be recreated
        timeit                      Boolean, if True, elapsed time will print
        save                        Boolean, if True, dataframe will save. Nice
                                    to set to False when testing.

    Output:
        df                          DataFrame, gets saved to pickle file.

    """
    if init:
        df = pd.DataFrame(columns=['Domain', 'Subdomain', 'Graph_order',
                'fp_gml', 'Weighted', 'Directed', 'Bipartite', 'Multigraph', 'Multiplex'])
    else:
        df = pd.read_pickle('newgmlcatalog.p')
    start = time.time()
    df = buildGMLcatalog(gmldirpath, df, avoiddirs, overwrite)
    elapsed = time.time()-start
    if timeit:
        print 'Elapsed time is %s seconds' %elapsed
    if save:
        df.to_pickle('newgmlcatalog.p')


if __name__ == '__main__':
    # path to gmls
    localpath = '/Users/annabroido/CU/gmls'
    externalpath = '/Volumes/External/gmls'
    avoiddirs = [ 'n7', 'n8', 'n9']
    # avoiddirs=[]
    #df = buildGMLcatalog(externalpath, df, avoiddirs, overwrite=True)
    preprocess(externalpath, timeit=True, avoiddirs=avoiddirs, save=True)
