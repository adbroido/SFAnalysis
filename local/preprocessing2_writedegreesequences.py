import numpy as np
import igraph
import pickle
import pandas as pd
import sortgmls as sg
import collections

import traceback
import warnings
import sys

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    traceback.print_stack()
    log = file if hasattr(file,'write') else sys.stderr
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

warnings.showwarning = warn_with_traceback

def writeerror(errormessage, size):
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

def readdeg(g, fp, degdir, analysis, namekey='', bipkey=0, weighkey=0, dirkey=0, mgkey=0, mpkey=0, compkey='entire'):
    """ Reads in an igraph object and writes the degree sequence to a text file.
    Assumes that g has been processed already and is not weighted, multiplex,
    bipartite, or directed.

    Input:
        g                       igraph object, represents one network
        fp                      path to GML file where a version of g lives
        degdir                  directory, where to store the degree sequences
        direction               in or out, indicates whether to take in or out degree

    Output:
        degfile                 output is written to a text file, not returned

    """
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
    if np.mean(deg) > np.sqrt(len(deg)):
        errormessage = "%s is too dense \n" %fn
        writeerror(errormessage, 'big')
    # check if mean degree is too small and don't write to file if so
    elif np.mean(deg) < 2:
        errormessage = "%s mean degree is too small \n" %fn
        writeerror(errormessage, 'small')
    else:
        count_dict = sorted(collections.Counter(deg).items())
        numedges = g.ecount()
        df = pd.DataFrame(count_dict, columns = ['xvalue', 'counts'])
        csvfile = degdir+fn
        df.to_csv(csvfile, index=False)
        # add new line if not already in index
        if fn not in analysis.index:
            analysis.loc[fn] = ''
        analysis.loc[fn]['Domain'] = domain
        analysis.loc[fn]['Subdomain'] = subdomain
        analysis.loc[fn]['fp_gml'] = fp
        analysis.loc[fn]['Graph_order'] = gsize
        analysis.loc[fn]['num_edges'] = numedges
        analysis.loc[fn]['Weighted'] = weighkey
        analysis.loc[fn]['Directed'] = dirkey
        analysis.loc[fn]['Bipartite'] = bipkey
        analysis.loc[fn]['Multigraph'] = mgkey
        analysis.loc[fn]['Multiplex'] = mpkey
        analysis.loc[fn]['Component'] = compkey

def checkconnected(g,fp,degdir, analysis, namekey='', bipkey=0, weighkey=0, dirkey=0, mpkey=0, mgkey=0):
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
    keyV = [('in', '_directedin'), ('out','_directedout'), ('total', '_directedtotal')]
    for dirkey, dirnamekey in keyV:
        newnamekey = namekey+dirnamekey
        checkconnected(g,fp,degdir, analysis, namekey=newnamekey, bipkey=bipkey, weighkey=weighkey, dirkey=dirkey, mgkey=0, mpkey=mpkey)

def find_threshold(weights, target_num_edges, left=0):
    uniqueweights = np.unique(weights)
    numbins = len(uniqueweights)
    done = False
    right = numbins-1
    while not done:
        mid = int(np.floor((left+right)/2))
        lenwtail = np.sum(weights>uniqueweights[mid])
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
    edgeseq = g.es(weight_gt=thresh)
    subg = g.subgraph_edges(edgeseq)
    if subg.is_directed():
        processdirected(subg,fp,degdir,analysis, namekey=namekey,weighkey=weighkey)
    else:
        checkconnected(subg,fp,degdir,analysis, namekey=namekey,weighkey=weighkey)

def processweighted(g, fp, degdir, analysis):
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
    if sg.weighted(g, fp)==1 or sg.multigraph(g)==1:
        processmultigraph(g, fp, degdir, analysis, namekey, mpkey=0, bipkey=bipkey)
    elif g.is_directed():
        processdirected(g,fp, degdir, analysis, namekey=namekey, mpkey=0, bipkey=bipkey)
    else:
        checkconnected(g,fp,degdir,analysis, namekey=namekey, mpkey=0, bipkey=bipkey)


def processbipartite(g, fp, degdir, analysis, namekey='', mpkey=0):
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

def processgraphs(catalog, degdir, analysis, overwrite=False):
    fpV = catalog['fp_gml']
    # only add new files
    if overwrite == False:
        fpV = set(fpV).difference(set(analysis['fp_gml']))
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
    # where to store deg seqs
    degdir = '/Users/annabroido/Dropbox/Research/LRTAnalysis/degreesequences/'
    catalog = pd.read_pickle('/Users/annabroido/Dropbox/Research/LRTAnalysis/LRTAnalysis/analysis/gmlcatalog.p')
    # analysis = pd.DataFrame(columns=['Domain', 'Subdomain', 'num_edges',
    #                                  'Graph_order', 'Weighted', 'Directed',
    #                                  'Bipartite', 'Multigraph', 'Multiplex',
    #                                  'Component', 'fp_gml', 'n', 'alpha',
    #                                  'xmin', 'ntail', 'Lpl', 'ppl', 'dexp',
    #                                  'dln', 'dstrexp', 'dplwc'])
    analysis = pd.read_pickle('/Users/annabroido/Dropbox/Research/LRTAnalysis/LRTAnalysis/analysis/analysis.p')
    # trim catalog to relevant entries
    catalog = catalog.query('Graph_order > 5 & Graph_order < 7')
    #fp = '/Volumes/Durodon/gmls/Biological/Food_web/n2/Aishihik_Lake_host-parasite_web_Aishihik_Lake_host-parasite_web_Biological_Food_web_n2.gml'
    processgraphs(catalog,degdir,analysis, overwrite=True)
    analysis.to_pickle('/Users/annabroido/Dropbox/Research/LRTAnalysis/LRTAnalysis/analysis/subanalysis_6.p')
