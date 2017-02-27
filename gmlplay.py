import igraph
import glob
import numpy as np

def myunion(graphs):
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




gmldpn2 = '/Users/annabroido/CU/testgmls/Social/Political/n2/'
gmldpn3 = '/Users/annabroido/CU/testgmls/Social/Political/n3/'


senatefpV = glob.glob(gmldpn2+'*enate*.gml')
senatefpV += glob.glob(gmldpn3+'*enate*.gml')
housefpV = glob.glob(gmldpn2+'*ouse*.gml')
housefpV += glob.glob(gmldpn3+'*ouse*.gml')

graphs = []
for fp in senatefpV:
    # read graph
    g = igraph.read(fp)
    # get snapshot id
    fn = fp.split("/")[-1]
    fnsplit = fn.split("_")
    ind = fnsplit.index("Senate")+1
    snap = fnsplit[ind]
    # add to edges
    g.es['snap'] = snap
    graphs.append(g)
