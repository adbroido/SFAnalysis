import fit
import lrt
import importfiles as im
import sortgmls as sg
import numpy as np
import pandas as pd
import igraph
import os
import collections



def buildGMLcatalog(gml_dir):
    """ Walks through the subdirectories of a root to find all gml files, then
    catalogs the relevant information about the contained networks.

    Input:
        gmldirpath              string, path to the root directory where gmls are



    Output:
        df                      DataFrame, catalog of the existing gml files

    """
    df = pd.DataFrame(columns=['fp_gml', 'Weighted', 'Directed', 'Bipartite',
                               'Multigraph', 'Multiplex'])
    # make list of file paths to gmls
    fpV = []
    for root, dirs, files in os.walk(gml_dir):
        for name in files:
            # leave out the bipartite projections so we can make our own
            if name.endswith('.gml'):
                fpV.append(os.path.join(root, name))
    # create the catalog
    for fp in fpV:
        g = igraph.read(fp)
        splitfp = fp.split('/')
        name = splitfp[-1]
        # add new row or overwrite existing row
        df.loc[name] = np.nan
        df.loc[name]['fp_gml'] = fp
        df.loc[name]['Weighted'] = sg.weighted(g, fp)
        df.loc[name]['Directed'] = sg.directed(g)
        df.loc[name]['Bipartite'] = sg.bipartite(g, fp)
        df.loc[name]['Multigraph'] = sg.multigraph(g)
        df.loc[name]['Multiplex'] = sg.multiplex(g)
        if (df.loc[name] == 'error').any():
            # this catches bad bipartite gmls
            df = df.drop(name)
            print('dropping {} from the considered gmls'.format(name))
    return df

"""
The functions here are designed to take cataloged gml files, extract
corresponding simple degree sequences, and store basic information about each
new simple graph and what path created it. The function processgraphs() is the
only one intended to be called directly. It requires as input a catalog of gmls
with file paths, and structural information (directed, weighted, multiplex,
multigraph, bipartite).
"""


def writeerror_deg(errormessage, size):
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

def readdeg(g, fp, degdir, analysis, namekey='', bipkey=0, weighkey=0, dirkey=0, mgkey=0, mpkey=0):
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
    if len(splitfp)>1:
        domain = splitfp[-4]
        subdomain = splitfp[-3]
        gsize = int(splitfp[-2][1:]) # cut out the 'n'
        gmlname = splitfp[-1]
    else:
        domain = 'na'
        subdomain = 'na'
        gsize = 'na'
        gmlname = fp
    fn = gmlname+namekey + 'distribution.txt'
    # check if degree sequence is too dense and don't write to file if so
    meandeg = np.mean(deg)
    if meandeg> np.sqrt(len(deg)):
        errormessage = "%s is too dense \n" %fn
        print(errormessage)
        writeerror_deg(errormessage, 'big')
    # check if mean degree is too small and don't write to file if so
    elif meandeg < 2:
        errormessage = "%s mean degree is too small \n" %fn
        print(errormessage)
        writeerror_deg(errormessage, 'small')
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

def processdirected(g, fp, degdir, analysis,  namekey='', bipkey=0, mpkey=0, weighkey=0, mgkey=0):
    """ Processes a directed graph. The graph is split into three: in-degree,
    out-degree, and total degree.

    """
    keyV = [('in', '_directedin'), ('out','_directedout'), ('total', '_directedtotal')]
    for dirkey, dirnamekey in keyV:
        newnamekey = namekey+dirnamekey
        readdeg(g,fp,degdir, analysis, namekey=newnamekey, bipkey=bipkey, weighkey=weighkey, dirkey=dirkey, mgkey=0, mpkey=mpkey)

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
        readdeg(subg,fp,degdir,analysis, namekey=namekey,weighkey=weighkey)

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
        readdeg(g,fp,degdir,analysis, namekey=namekey, mpkey=mpkey, bipkey=bipkey, mgkey=mgkey, weighkey=weighkey)



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
        readdeg(g,fp,degdir,analysis, namekey=namekey, mpkey=0, bipkey=bipkey)


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
    print fp
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
                readdeg(graph, fp,degdir,analysis, namekey=namekey, mpkey=mpkey)
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
                readdeg(graph, fp,degdir,analysis, namekey=namekey, mpkey=mpkey)

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
            readdeg(graph, fp, degdir, analysis, namekey=namekey, mpkey=mpkey)


def write_degree_sequences(gml_dir, deg_dir):
    gml_df = buildGMLcatalog(gml_dir)
    fpV = gml_df['fp_gml']
    analysis_df = pd.DataFrame(columns=['num_edges', 'Weighted', 'Directed',
                                         'Bipartite', 'Multigraph', 'Multiplex',
                                         'fp_gml', 'n', 'alpha', 'xmin','ntail',
                                         'Lpl', 'ppl', 'dexp', 'dln', 'dstrexp',
                                         'dplwc', 'meandeg'])
    for fp in fpV:
        g = igraph.read(fp)
        #### find what kind of graph this is (follow hierarchical ordering of types)
        # check first for multiplex
        row = gml_df[gml_df.fp_gml==fp]
        if row['Multiplex'].item() == 1:
            processmultiplex(g,fp, deg_dir, analysis_df)
        elif row['Bipartite'].item() == 1:
            processbipartite(g, fp, deg_dir,analysis_df)
        elif row['Multigraph'].item() == 1:
            processmultigraph(g, fp, deg_dir,analysis_df)
        elif row['Weighted'].item() == 1:
            processweighted(g, fp, deg_dir,analysis_df)
        elif row['Directed'].item() == 1:
            processdirected(g, fp, deg_dir,analysis_df)
        else:
            readdeg(g, fp,deg_dir,analysis_df)
    return analysis_df

def organize_degree_sequences(deg_dir):
    fnV = [file for file in os.listdir(deg_dir) if file.split('.')[-1] in
                                                                ['txt', 'csv']]
    analysis_df = pd.DataFrame(columns=['fp_gml', 'n', 'alpha', 'xmin','ntail',
                                        'Lpl', 'ppl', 'dexp', 'dln', 'dstrexp',
                                         'dplwc', 'meandeg'], index = fnV)

    for fn in fnV:
        analysis_df.loc[fn]['fp_gml'] = 'na'
    return analysis_df

def writeerror_analysis(errormessage):
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

def writeerror_lrt(errormessage):
    errorfp = 'lrterror.txt'
    known = False
    f = open(errorfp, 'a+')
    f.seek(0)
    for line in f:
        if line == errormessage:
            known = True
    if not known:
        f.write(errormessage)
    f.close()


def analyze_degree_sequences(deg_dir, analysis):
    for fn in analysis.index:
        fp = deg_dir + fn
        x = im.readdata(fp)
        # note if there is a problem with the file
        n = len(x)
        if np.mean(x) < 2 or np.mean(x) > np.sqrt(n):
            errormessage = "%s has a bad mean degree \n" %fp
            writeerror_analysis(errormessage)
        # catch for trivial degree sequences
        elif len(np.unique(x)) == 1:
            errormessage = "%s contains only one unique value \n" %fp
            writeerror_analysis(errormessage)
        else:
            if analysis.loc[fn]['ppl']=='' or overwrite == True:
                analysis.loc[fn]['n'] = n
                [alpha, xmin, ntail,  L, ks] = fit.pl(x)
                p = fit.plpval(x,alpha, xmin, ks)
                analysis.loc[fn]['alpha'] = alpha
                analysis.loc[fn]['xmin'] = xmin
                analysis.loc[fn]['ntail'] = ntail
                analysis.loc[fn]['Lpl'] = L
                analysis.loc[fn]['ppl'] = p
            if analysis.loc[fn]['dexp']=='' or overwrite == True:
                # compare the alternative distributions
                xmin = analysis.loc[fn]['xmin']
                alpha = analysis.loc[fn]['alpha']
                x = x[x>=xmin]
                # compare the non-nested alternatives, return the decisions for each
                decisionthresh = 0.1
                [dexp, dln, dstrexp] = lrt.nonnested(x,alpha, decisionthresh)
                if dexp == 2:
                    errormessage = "Exponential didn't converge for %s \n" %fp
                    writeerror_lrt(errormessage)
                if dln == 2:
                    errormessage = "Log-normal didn't converge for %s \n" %fp
                    writeerror_lrt(errormessage)
                if dstrexp == 2:
                    errormessage = "Stretched exponential didn't converge for %s \n" %fp
                    writeerror_lrt(errormessage)
                # fit the nested alternatives
                dplwc = lrt.nested(x, alpha, decisionthresh)
                if dplwc == 2:
                    errormessage = "PLWC didn't converge for %s \n" %fp
                    writeerror_lrt(errormessage)
                # update dataframe
                analysis.loc[fn]['dexp'] = dexp
                analysis.loc[fn]['dln'] = dln
                analysis.loc[fn]['dstrexp'] = dstrexp
                analysis.loc[fn]['dplwc'] = dplwc
    return analysis

""" Helper functions for categorizing networks into scale-free types"""
def test_strong(rows):
    S1 = False # strongest
    S2 = False # strong
    SA = False
    n = len(rows)
    strong = 0
    strong_alone = 0
    for ind, row in rows.iterrows():
        if row.ppl>0.1 and row.ntail >= 50 and row.alpha < 3 and row.alpha > 2:
            strong_alone += 1
            if row.dexp >-1 and row.dln>-1  and row.dstrexp >-1  and row.dplwc >-1:
                strong += 1
    if strong_alone >= 9.*n/10:
        SA = True
    if strong >= n/2.:
        S2 = True
    if SA == True and strong >= 95.*n/100:
        S1 = True
    return (S1, S2)

def test_weak(rows):
    W = False
    West = False
    SW = False
    n = len(rows)
    weak = 0
    weakest = 0
    sweak = 0
    for ind, row in rows.iterrows():
        if row.ppl>0.1:
            weakest += 1
            if row.ntail>=50:
                weak += 1
        if row.dexp >-1 and row.dln>-1  and row.dstrexp >-1  and row.dplwc >-1:
            sweak += 1
    if weak >= n/2.:
        W = True
    if weakest >= n/2.:
        West = True
    if sweak >= n/2.:
        SW = True
    return (W, West, SW)

def test_strong_any(rows):
    ''' just need a minimum of one dataset'''
    S = False # strongest/strong
    n = len(rows)
    strong_alone = 0 #noplwc
    strong = 0
    for ind, row in rows.iterrows():
        if row.ppl>0.1 and row.ntail >= 50 and row.alpha <3 and row.alpha > 2:
            if row.dexp >-1 and row.dln>-1  and row.dstrexp >-1 and row.dplwc >-1:
                strong += 1
    if strong>= 1:
        S = True
    return (S)

def test_weak_any(rows):
    ''' just need a minimum of one dataset'''
    W = False
    West = False
    SW = False
    n = len(rows)
    weak = 0
    weakest = 0
    sweak = 0
    for ind, row in rows.iterrows():
        if row.ppl>0.1:
            weakest += 1
            if row.ntail>=50:
                weak += 1
        if row.dexp >-1 and row.dln>-1  and row.dstrexp >-1  and row.dplwc >-1:
            sweak += 1
    if weak >= 1:
        W = True
    if weakest >= 1:
        West = True
    if sweak >= 1:
        SW = True
    return (W, West, SW)


def categorize_networks(df, permissive=False):
    unique_datasets = np.unique(df.fp_gml)
    if permissive:
        hyps = pd.DataFrame(columns = ['Strongest', "Strong", "Weak", "Weakest",
                                       "Super_Weak", "Weak_Any", "Weakest_Any",
                                       "Super_Weak_Any", "median_alpha", "n",
                                       "median_ntail"], index=unique_datasets )
    else:
        hyps = pd.DataFrame(columns = ['Strongest', "Strong", "Weak", "Weakest",
                                       "Super_Weak", "median_alpha", "n",
                                       "median_ntail"], index=unique_datasets )
    for i, dataset in enumerate(unique_datasets):
        query = "fp_gml == '%s'" %dataset
        rows = df.query(query)
        [S1, S2] = test_strong(rows)
        hyps.loc[dataset]['Strongest'] = S1
        hyps.loc[dataset]['Strong'] = S2
        [weak,weakest,superweak] = test_weak(rows)
        hyps.loc[dataset]['Weak'] = weak
        hyps.loc[dataset]['Weakest'] = weakest
        hyps.loc[dataset]['Super_Weak'] = superweak
        if permissive:
            S = test_strong_any(rows)
            hyps.loc[dataset]['Strong_Any'] = S
            [W, West, SW] = test_weak_any(rows)
            hyps.loc[dataset]['Weak_Any'] = W
            hyps.loc[dataset]['Weakest_Any'] = West
            hyps.loc[dataset]['Super_Weak_Any'] = SW
        hyps.loc[dataset]['median_alpha'] = np.median(rows.alpha)
        hyps.loc[dataset]['n'] = np.max(rows.n)
        hyps.loc[dataset]['median_ntail'] = np.median(rows.ntail)
    return hyps
