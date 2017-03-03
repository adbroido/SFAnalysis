import igraph
import glob
import numpy as np

def union(graphs):
    if not (all([graph.is_directed() for graph in graphs]) or all([not graph.is_directed() for graph in graphs])):
        print "You can't mix directed and undirected graphs"
    v = np.max([graph.vcount() for graph in graphs])
    g = igraph.Graph(n=v, directed=graphs[1].is_directed())
    for graph in graphs:
        attrnames = graph.es.attributes()
        start = g.ecount()
        g.add_edges([e.tuple for e in graph.es])
        stop = g.ecount()
        for att in attrnames:
            g.es[start:stop][att] = graph.es[att]
    return g

def combine_gmls(fpV, splitkey, typekey='snap', comment=''):
    graphs = []
    for fp in fpV:
        # read graph
        g = igraph.read(fp)
        # get snapshot id
        fn = fp.split("/")[-1]
        fnsplit = fn.split("_")
        ind = fnsplit.index(splitkey)+1
        snap = fnsplit[ind]
        # add to edges
        g.es[typekey] = snap
        graphs.append(g)
    # combine graphs
    g = union(graphs)
    # write to gml file
    nameV = fn.split(snap)
    gmlname = nameV[0] + 'Union' + nameV[1]
    comment = 'Union of the ' + nameV[0]+snap+' files. ' + comment
    g.write_gml(gmlname, comment)

def congress():
    """
    Union the Congress graphs - One for Senate and one for House
    """
    # file paths to the gml directories
    gmldpn2 = '/Users/annabroido/CU/testgmls/Social/Political/n2/'
    gmldpn3 = '/Users/annabroido/CU/testgmls/Social/Political/n3/'

    # files paths, separated by type
    senatefpV = glob.glob(gmldpn2+'*enate*.gml')
    senatefpV += glob.glob(gmldpn3+'*enate*.gml')
    housefpV = glob.glob(gmldpn2+'*ouse*.gml')
    housefpV += glob.glob(gmldpn3+'*ouse*.gml')

    combine_gmls(senatefpV, splitkey='Senate', typekey='snap')
    combine_gmls(housefpV, splitkey='House', typekey='snap')

def norwegian():
    """
    Union the Congress graphs - One for Senate and one for House
    """
    # file paths to the gml directories
    gmldp = '/Users/annabroido/CU/testgmls/Social/Affiliation/n4/'
    # files paths, separated by type
    fpV = glob.glob(gmldp+'*net2mode*.gml')

    combine_gmls(fpV, splitkey='net2mode', typekey='snap')

#congress()
norwegian()
