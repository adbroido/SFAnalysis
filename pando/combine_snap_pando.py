import os
import glob
import igraph
import numpy as np
import itertools

root_directory = "/Users/anbr3575/LRTAnalysis/gmls/Technological/Communication/"
text_to_match = "Route_Views_AS*1999-2000*"
splitkey = "1999-2000"

paths = []
for nodedir in ["n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8"]:
    path_to_search = root_directory + nodedir + "/"
    if os.path.isdir(path_to_search):
        local_paths = glob.glob(path_to_search + text_to_match + "*.gml")
        if local_paths:
            paths += local_paths
            order = nodedir

full_g = igraph.read(paths[0])
snap = paths[0].split('/')[-1].split(splitkey)[1].split('_')[1]
full_g.es['snap'] = snap

for partial_ind in range(1,len(paths)):
    print partial_ind
    partial_g = igraph.read(paths[partial_ind])

    full_nodes = [node.attributes() for node in full_g.vs]
    partial_nodes = [node.attributes() for node in partial_g.vs]

    #find new nodes
    new_nodes = list(itertools.ifilterfalse(lambda x: x in full_nodes, partial_nodes))

    # add new nodes to full_g
    for node in new_nodes:
        new_node_ind = len(full_g.vs)
        full_g.add_vertices(1)
        full_g.vs[new_node_ind]['id'] = node['id']
    full_nodes = [node.attributes() for node in full_g.vs]

    # get snapshot id
    snap = paths[partial_ind].split('/')[-1].split(splitkey)[1].split('_')[1]

    # add edges
    partial_g.es['snap'] = snap
    start = full_g.ecount()
    attrnames = full_g.es.attributes()
    for edge_ind, edge  in enumerate(partial_g.es()):
        (local_source, local_target) = edge.tuple
        full_source = full_nodes.index(partial_nodes[local_source])
        full_target = full_nodes.index(partial_nodes[local_target])
        full_g.add_edges([(full_source,full_target)])
        for att in attrnames:
            full_g.es[start+edge_ind][att] = edge[att]

# # write to gml file
fnsplit = fn.split("_" + snap + "_")
orderofmag = int(np.ceil(np.log10(full_g.vcount())))
order = "n%s" %orderofmag
out_gml_title = "/Users/anbr3575/" + fnsplit[0] + "_Union_" + fnsplit[1].split(".")[0][:-2] + order + ".gml"
print out_gml_title
comment = "Union of the " + fnsplit[0][:-1] + " files."
full_g.write_gml(out_gml_title, comment)
