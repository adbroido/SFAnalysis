import pandas as pd
import numpy as np
import itertools as it

def readdataold(fp):
    """ Reads in a datafile. Assumes the first two rows are headers we may wish
    to keep, and that the remaining rows are integers. Uses generators to
    efficiently loop through the rows and avoids intermediate step of saving the
    entire list of integers as strings (roughly 5x more memory). About a 5x
    improvement in performance over np's genfromtxt, also saves the first two
    rows as strings.

    Input:
        fp                     string, filepath (must be in python's path).
                               Assumes the first 2 rows are single-word strings
                               we want to keep, the third row is an integer we
                               want to keep, and that all subsequent rows are
                               data (integers).

    Output:
        domain                 string, first row, corresponds to domain of graph
        subdomain              string, second row, subdomain of the graph
        edges                  int, third row, number of edges in graph
        gsize                  string, order of number of nodes (ie n4 = 10^4)
        data                   ndarray, ndim = 1, dtype = integer.
    """

    # create a generator that strips all the extra whitespace
    lines = (line.rstrip() for line in open(fp))
    # the first four rows are headers we want to keep
    domain, subdomain, edgesstr, gsize = [text for text in it.islice(lines, 4)]
    # convert the third one to an integer
    edges = int(edgesstr)
    # the rest of the data are integers that form the dataset
    data = np.array([int(row) for row in lines])
    return [domain, subdomain, edges, gsize, data]

def readdata(fp):
    """ Reads in a datafile.

    Input:
        fp                     string, filepath (must be in python's path).

    Output:
        data                   ndarray, ndim = 1, dtype = integer.
    """

    df = pd.read_csv(fp)
    data = np.concatenate([np.array([df.xvalue[i] for repeat in range(df.counts[i])]) for i in range(len(df))])
    return data

if __name__ == '__main__':
    degseqdirpath = '/Users/annabroido/Dropbox/Research/LRTAnalysis/Code/degreesequences/*done*/'
    samplefile = '9-11_terrorist_network_hijacker_associations_Social_Offline_n2.gmldistribution.txt'
    fp = degseqdirpath + samplefile
    h1, h2, e, data = readdata(fp)
