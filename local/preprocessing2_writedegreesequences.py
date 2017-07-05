import numpy as np
import igraph
import pickle
import pandas as pd
import sortgmls as sg
import collections

"""
The functions here are designed to take cataloged gml files, extract
corresponding simple degree sequences, and store basic information about each
new simple graph and what path created it. The function processgraphs() is the
only one intended to be called directly. It requires as input a catalog of gmls
with file paths, and structural information (directed, weighted, multiplex,
multigraph, bipartite).
"""


def writeerror(errormessage, size):
    """ Write to error files for degree sequences whose mean is either too big
    or too small to be considered.

    Input:
        errormessage            string, gets written as a single line in the
                                error file.
        size                    string, either 'small' or 'big'. Designates
                                which error file to write to.
    """
    if size=='small':
        errorfp = 'degseqerror_small.txt'
    elif size=='big':
        errorfp = 'degseqerror_big.txt'
    known = False
    f = open(errorfp, 'a+')
    f.seek(0)
    for line in f:
        if line == errormessage:
            known = True
    if not known:
        f.write(errormessage)
    f.close()

def readdeg(g, fp, degdir, analysis, namekey='', bipkey=0, weighkey=0, dirkey=0, mgkey=0, mpkey=0, compkey):
    """ Reads in an igraph object and writes the degree sequence to a text file.
    Assumes that g has been processed already and is simple or directed only.
    A new row with basic information about the graph is added to analysis.

    Input:
        g                       igraph object, represents one graph
        fp                      path to GML file where a version of g lives
        degdir                  directory, where to store the degree sequences
        analysis                DataFrame, table of results indexed by degree seq
        namekey                 string, gets added to the degree sequence
                                filename to indicate the path that created the
                                sequence from the original network
        bipkey                  int or string, indicates type of bipartite
                                splitting. 0 for not bipartite.
        weighkey                int or string, indicates weighting threshold or
                                simplification algorithm. 0 for not weighted.
        dirkey                  int or string, indicates direction of edges for
                                splitting. 0 for not directed.
        mgkey                   int or string, indicates simplified from
                                multigraph 0 for not multigraph.
        mpkey                   int or string, indicates edge type that generated
                                this subgraph. 0 for not weighted.
        compkey                 string, indicates weighting threshold or
                                simplification algorithm. 0 for not weighted.


    Output:
        degfile                 output is written to a text file, not returned

    """
    # split by degree
    if dirkey == 0 or dirkey == 'total':
        deg = g.degree()
    elif dirkey == 'in':
        deg = g.indegree()
    elif dirkey == 'out':
        deg = g.outdegree()
    else:
        print 'something is wrong with your dirkey'
    # get file name
    splitfp = fp.split('/')
    domain = splitfp[-4]
    subdomain = splitfp[-3]
    gsize = int(splitfp[-2][1:]) # cut out the 'n'
    gmlname = splitfp[-1]
    fn = gmlname+namekey + 'distribution.txt'
    # check if degree sequence is too dense and don't write to file if so
    meandeg = np.mean(deg)
    if meandeg> np.sqrt(len(deg)):
        errormessage = "%s is too dense \n" %fn
        writeerror(errormessage, 'big')
    # check if mean degree is too small and don't write to file if so
    elif meandeg < 2:
        errormessage = "%s mean degree is too small \n" %fn
        writeerror(errormessage, 'small')
    else:
        # write degree sequence file. Each row is xvalue,count
        count_dict = sorted(collections.Counter(deg).items())
        numedges = g.ecount()
        df = pd.DataFrame(count_dict, columns = ['xvalue', 'counts'])
        csvfile = degdir+fn
        df.to_csv(csvfile, index=False)
        # add new line to the data frame if not already in index
        if fn not in analysis.index:
            analysis.loc[fn] = ''
        analysis.loc[fn]['Domain'] = domain
        analysis.loc[fn]['Subdomain'] = subdomain
        analysis.loc[fn]['fp_gml'] = fp
        analysis.loc[fn]['Graph_order'] = gsize
        analysis.loc[fn]['num_edges'] = numedges
        analysis.loc[fn]['meandeg'] = meandeg
        analysis.loc[fn]['Weighted'] = weighkey
        analysis.loc[fn]['Directed'] = dirkey
        analysis.loc[fn]['Bipartite'] = bipkey
        analysis.loc[fn]['Multigraph'] = mgkey
        analysis.loc[fn]['Multiplex'] = mpkey
        analysis.loc[fn]['Component'] = compkey

def checkconnected(g,fp,degdir, analysis, namekey='', bipkey=0, weighkey=0, dirkey=0, mpkey=0, mgkey=0):
    """ Checks if a graph is connected. If not, pull out 2 degree sequences:
    entire graph and jsut the largest component.

    """
    # if g is not connected, pull out the largest component
    if g.is_connected():
        compkey = 'connected'
        readdeg(g,fp,degdir, analysis, namekey=namekey, bipkey=bipkey, weighkey=weighkey, dirkey=dirkey, mpkey=mpkey, mgkey=mgkey, compkey=compkey)
    else:
        compkey = 'entire'
        readdeg(g,fp,degdir, analysis, namekey=namekey, bipkey=bipkey, weighkey=weighkey, dirkey=dirkey, mpkey=mpkey, mgkey=mgkey, compkey=compkey)
        clust = g.clusters() # note that the components function seems identical
        largecomp = g.subgraph(clust[0])
        newnamekey = namekey+ '_largestcomp'
        compkey = 'largest'
        readdeg(largecomp,fp,degdir, analysis, namekey=newnamekey, bipkey=bipkey, weighkey=weighkey, dirkey=dirkey, mpkey=mpkey, mgkey=mgkey, compkey=compkey)

def processdirected(g, fp, degdir, analysis,  namekey='', bipkey=0, mpkey=0, weighkey=0, mgkey=0):
    """ Processes a directed graph. The graph is split into three: in-degree,
    out-degree, and total degree.

    """
    keyV = [('in', '_directedin'), ('out','_directedout'), ('total', '_directedtotal')]
    for dirkey, dirnamekey in keyV:
        newnamekey = namekey+dirnamekey
        checkconnected(g,fp,degdir, analysis, namekey=newnamekey, bipkey=bipkey, weighkey=weighkey, dirkey=dirkey, mgkey=0, mpkey=mpkey)

def find_threshold(weights, target_num_edges, left=0):
    """ Given a target number of edges in a new subgraph and a list of current
    edge weights, finds the threshold needed to achieve this target edge count.
    Uses bisection to find the threshold.

    Input:
        weights             ndarray, edge weights in the order corresponding to
                            the order of edges
        target_num_edges    float, threshold for keeping weights. Edges with
                            weight < thresh are discarded
        left                int, index of smallest weight to consider. Included
                            as an argument to speed up the search when we know
                            weaker thresholds.

    Output:
        thresh              float, threshold weight
        mid                 int, index of threshold weight

    """
    # ordered list of weights without duplicates
    uniqueweights = np.unique(weights)
    # we will bin weights into unique groups
    numbins = len(uniqueweights)
    done = False
    right = numbins-1
    while not done:
        # find midpoint
        mid = int(np.floor((left+right)/2))
        # number of weights above the midpoint
        lenwtail = np.sum(weights>uniqueweights[mid])
        # number of unique weights above mid
        lenutail = np.sum(uniqueweights>uniqueweights[mid])
        if lenwtail > target_num_edges:
            left = mid+ 1
        elif lenwtail < target_num_edges:
            right = mid- 1
        if left>right or lenwtail==target_num_edges or lenutail==1:
            done = True
    thresh= uniqueweights[mid]
    return thresh, mid

def oneweighted(g, fp, degdir, analysis, weights, thresh, namekey, weighkey):
    """ Processes a single weighted graph. Pulls out the list of edges with
    weight above the minimum threshold "thresh", creates a new subgraph from
    these edges, then checks for directedness and connectedness.

    Input:
        weights             ndarray, edge weights in the order corresponding to
                            the order of edges
        thresh              float, threshold for keeping weights. Edges with
                            weight < thresh are discarded
        weighkey            string, indicates which thresholding algorithm was
                            used

    """
    edgeseq = g.es(weight_gt=thresh)
    subg = g.subgraph_edges(edgeseq)
    if subg.is_directed():
        processdirected(subg,fp,degdir,analysis, namekey=namekey,weighkey=weighkey)
    else:
        checkconnected(subg,fp,degdir,analysis, namekey=namekey,weighkey=weighkey)

def processweighted(g, fp, degdir, analysis):
    """ Processes a weighted graph. This is only for graphs that are not
    multigraph, bipartite, or multiplex. The graph is split into three by
    thresholding on weight in different ways.

    """
    if not g.is_weighted():
        g.es['weight'] = g.es['value']
    weights = np.asarray(g.es['weight'])
    n = g.vcount()

    # w1: want <k>=sqrt(n), so m=(1/2)n^(3/2)
    namekey = '_weighted1'
    weighkey = 'w1'
    target_num_edges = (float(n)**(1.5))/2
    thresh1, ind1 = find_threshold(weights, target_num_edges)
    oneweighted(g, fp, degdir, analysis, weights,thresh1, namekey=namekey, weighkey=weighkey)

    # w2: want <k> in between 2 and sqrt(n)
    namekey = '_weighted2'
    weighkey = 'w2'
    target_num_edges = (float(n)**(float(5)/4))/2
    thresh2, ind2 = find_threshold(weights, target_num_edges, left=ind1)
    # if the threshold is the same, don't rewrite the deg seqs
    if ind2 > ind1:
        oneweighted(g, fp, degdir, analysis, weights,thresh2, namekey, weighkey)

    # w3: want <k>=2, so m=n
    namekey = '_weighted3'
    weighkey = 'w3'
    target_num_edges = n
    thresh3, ind3 = find_threshold(weights, target_num_edges, left=ind2)
    # if the threshold is the same, don't rewrite the deg seqs
    if ind3 > ind2:
        oneweighted(g, fp, degdir, analysis, weights,thresh3, namekey, weighkey)

def processmultigraph(g, fp, degdir, analysis, namekey='', mpkey=0, bipkey=0):
    """ Processes a multigraph or weighted graph by ignoring multiedges and
    weights, then sends the simplified graph on.

    """
    mgkey=0
    weighkey=0
    if sg.multigraph(g)==1:
        namekey += '_multigraphsimplified'
        mgkey = 'simplified'
    if sg.weighted(g)==1:
        namekey += '_weightedsimplified'
        weighkey = 'simplified'
    g.simplify()
    if sg.directed(g):
        processdirected(g,fp,degdir,analysis, namekey=namekey, mpkey=mpkey, bipkey=bipkey, mgkey=mgkey, weighkey=weighkey)
    else:
        checkconnected(g,fp,degdir,analysis, namekey=namekey, mpkey=mpkey, bipkey=bipkey, mgkey=mgkey, weighkey=weighkey)



def onebipartite(g,fp,degdir, analysis, namekey,mpkey, bipkey):
    """ Processes a single bipartite projection, sending it to the next step in
    the hierarchy towards writing degree seqeunces.

    Input:
        bipkey                  string, indicates the projection (a or b) that
                                generated this subgraph of the original.

    """
    if sg.weighted(g, fp)==1 or sg.multigraph(g)==1:
        processmultigraph(g, fp, degdir, analysis, namekey, mpkey=0, bipkey=bipkey)
    elif g.is_directed():
        processdirected(g,fp, degdir, analysis, namekey=namekey, mpkey=0, bipkey=bipkey)
    else:
        checkconnected(g,fp,degdir,analysis, namekey=namekey, mpkey=0, bipkey=bipkey)


def processbipartite(g, fp, degdir, analysis, namekey='', mpkey=0):
    """ Processes a bipartite graph. Splits into three graphs: a- and b-mode
    projections, and full graph. Sends each on along to the next step in the
    hierarchy towards writing degree seqeunces.

    Input:
        namekey                 string, gets updated at every step in the
                                hierarchy and in the end is added to the name of
                                the degree sequence file to keep track of what
                                path generated each sequence
        mpkey                   int or string, indicates the edge type that
                                generated this subgraph of the original. 0 if
                                the original network was not multiplex

    """
    # get node types and recast as integer types
    types = np.asarray(g.vs['type'])
    uniquetypes = np.unique(types)
    types[types==uniquetypes[0]] = 0
    types[types==uniquetypes[1]] = 1
    types = np.asarray([int(t) for t in types])
    g.vs['type'] = types

    a,b = g.bipartite_projection()

    # deal with a
    newnamekey = namekey+'_bipartitea'
    bipkey = 'a'
    onebipartite(a,fp,degdir, analysis, namekey=newnamekey,mpkey=mpkey,bipkey=bipkey)

    # deal with b
    newnamekey = namekey+'_bipartiteb'
    bipkey = 'b'
    onebipartite(b,fp,degdir,analysis, namekey=newnamekey,mpkey=mpkey,bipkey=bipkey)

    # deal with entire graph
    newnamekey = namekey+'_bipartitefull'
    bipkey = 'full'
    onebipartite(g,fp,degdir, analysis, namekey=newnamekey,mpkey=mpkey,bipkey=bipkey)

def processmultiplex(g, fp, degdir, analysis):
    """ Processes a multiplex graph. Splits along the edge types, so that each
    edge type gets its own new graph. Then sends each of these new graphs
    through the structural hierarchy and sends them along the appropriate path
    for further processing.

    Input:
        g                     igraph Graph object, known to be multiplex
        fp                    file path, leads to gml file

    """
    # project onto layers
    # pull out list of attributes
    attributes = g.es.attributes()
    # if there are multiple edge types, split on these.
    if len(attributes)>1:
        for att in attributes:
            # assume that these attribute values are weights, so should be numeric
            attkeyword = att+'_notin'
            # list of non-numeric values to avoid
            notthese = ['','Nan', 'n']
            # pull out edges that correspond to non empty weights
            edgeseq = g.es(**{attkeyword:notthese})
            # project onto the subgraph
            graph = g.subgraph_edges(edgeseq)
            namekey = '_multiplex'+att
            mpkey = 'sub_'+att
            if sg.bipartite(graph, fp)==1:
                processbipartite(graph, fp, degdir, analysis, namekey=namekey, mpkey=mpkey)
            elif sg.multigraph(graph)==1 or sg.weighted(graph)==1:
                processmultigraph(graph, fp, degdir, analysis, namekey=namekey, mpkey=mpkey)
            elif sg.directed(graph)==1:
                processdirected(graph, fp, degdir, analysis, namekey=namekey, mpkey=mpkey)
            else:
                checkconnected(g, fp,degdir,analysis, namekey=namekey, mpkey=mpkey)
    # If, however, there is one edge type, assume the split is in this, and
    # look at the values of this attribute
    else:
        att = attributes[0]
        # get just the unique values the attribute takes
        types = np.unique(g.es[att])
        # initilize empty list of edge seqs for each projection subgraph
        edgeseqs = []
        attkeyword = att+'_eq'
        for typ in types:
            edgeseqs.append((g.es(**{attkeyword:typ})))
        # project onto the subgraphs
        subgraphs = [g.subgraph_edges(edgeseq) for edgeseq in edgeseqs]
        # process all the subgraphs
        for i in range(len(subgraphs)):
            graph = subgraphs[i]
            namekey = '_multiplex'+str(types[i])
            mpkey = 'sub_'+str(i)
            if sg.bipartite(graph, fp)==1:
                processbipartite(graph, fp, degdir, analysis, namekey=namekey, mpkey=mpkey)
            elif sg.multigraph(graph)==1 or sg.weighted(graph)==1:
                processmultigraph(graph, fp, degdir, analysis, namekey=namekey, mpkey=mpkey)
            elif sg.directed(graph)==1:
                processdirected(graph, fp, degdir, analysis, namekey=namekey, mpkey=mpkey)
            else:
                checkconnected(g, fp,degdir,analysis, namekey=namekey, mpkey=mpkey)

        # process the union graph
        graph = g
        namekey = '_multiplexunion'
        mpkey = 'union'
        if sg.bipartite(graph, fp)==1:
            processbipartite(graph, fp, degdir, analysis, namekey=namekey, mpkey=mpkey)
        elif sg.multigraph(graph)==1:
            processmultigraph(graph, fp, degdir, analysis, namekey=namekey, mpkey=mpkey)
        elif sg.weighted(graph)==1:
            processweighted(graph, fp, degdir, analysis, namekey=namekey, mpkey=mpkey)
        elif sg.directed(graph)==1:
            processdirected(graph, fp, degdir, analysis, namekey=namekey, mpkey=mpkey)
        else:
            checkconnected(graph, fp, degdir, analysis, namekey=namekey, mpkey=mpkey)

def processgraphs(catalog, degdir, overwrite=False):
    """ Takes catalog of gml files and their structural properties(weighted,
    directed, bipartite, multigraph, multiplex) and sets them on the path to
    extract all relevant simple degree sequences from each.

    Input:
        catalog                     pandas DataFrame, indices are gml file names
        degdir                      file path, where degree seqeunce files will
                                    be stored
        overwrite                   boolean, if false only new files will be
                                    added to the analysis DataFrame. Otherwise
                                    the dataframe will be recreated from scratch.
                                    Default is False.

    """
    fpV = catalog['fp_gml']
    # only add new files
    if overwrite == False:
        analysis = pd.read_pickle('/Users/annabroido/Dropbox/Research/LRTAnalysis/LRTAnalysis/analysis/analysis.p')
        fpV = set(fpV).difference(set(analysis['fp_gml']))
    else:
        analysis = pd.DataFrame(columns=['Domain', 'Subdomain', 'num_edges',
                                         'Graph_order', 'Weighted', u'Directed',
                                         'Bipartite', 'Multigraph', 'Multiplex',
                                         'Component', 'fp_gml', 'n', 'alpha',
                                         'xmin', 'ntail', 'Lpl', 'ppl', 'dexp',
                                         'dln', 'dstrexp', 'dplwc', 'meandeg'])
    for fp in fpV:
        g = igraph.read(fp)
        #### find what kind of graph this is (follow hierarchical ordering of types)
        # check first for multiplex
        row = catalog[catalog.fp_gml==fp]
        if row['Multiplex'].item() == 1:
            processmultiplex(g,fp, degdir, analysis)
        elif row['Bipartite'].item() == 1:
            processbipartite(g, fp, degdir,analysis)
        elif row['Multigraph'].item() == 1:
            processmultigraph(g, fp, degdir,analysis)
        elif row['Weighted'].item() == 1:
            processweighted(g, fp, degdir,analysis)
        elif row['Directed'].item() == 1:
            processdirected(g, fp, degdir,analysis)
        else:
            checkconnected(g, fp,degdir,analysis)
    print 'All done!'



if __name__ == '__main__':
    # file path to degree sequences
    degdir = '/Users/annabroido/Dropbox/Research/LRTAnalysis/degreesequences/'
    # read in gml catalog to pandas DataFrame
    catalog = pd.read_pickle('/Users/annabroido/Dropbox/Research/LRTAnalysis/LRTAnalysis/analysis/gmlcatalog.p')
    # trim catalog to subset of entries. This step is optional
    catalog = catalog.query('Graph_order ==6')
    # run!
    processgraphs(catalog,degdir,analysis, overwrite=False)
    # save to pickle file
    analysis.to_pickle('/Users/annabroido/Dropbox/Research/LRTAnalysis/LRTAnalysis/analysis/analysis.p')
