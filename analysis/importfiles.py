import pandas as pd
import numpy as np


def readdata(fp):
    """ Reads in a datafile.

    Input:
        fp                      string, filepath to csv file (degree sequence).

    Output:
        data                    ndarray, ndim = 1, dtype = integer. Repeats each
                                xval as many times as indicated by counts
    """

    df = pd.read_csv(fp)
    data = np.concatenate([np.array([df.xvalue[i] for repeat in range(df.counts[i])]) for i in range(len(df))])
    return data
