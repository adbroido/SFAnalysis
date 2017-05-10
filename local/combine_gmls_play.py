import os
import glob
import igraph
import numpy as np

root_directory = "/Volumes/Durodon/clean_gmls/Social/Political/"
text_to_match = "*enate"
splitkey = "Senate"

paths = []
for nodedir in ["n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8"]:
    path_to_search = root_directory + nodedir + "/"
    if os.path.isdir(path_to_search):
        local_paths = glob.glob(path_to_search + text_to_match + "*.gml")
        if local_paths:
            paths += local_paths
            order = nodedir

graph0 = igraph.read(paths[0])
full_g = igraph.Graph(directed=graph0.is_directed())
full_node_dict = {}

# for fp in paths:
fp = paths[0]
g = igraph.read(fp)

#    DEAL WITH NODES

# make dictionary of indicies for nodes in g
local_node_dict = {}
g_vs = g.vs()
for node_ind, node in enumerate(g_vs):
    local_node_dict.update({node_ind:node.attributes()})
# add any new nodes to the full dictonary
for node in local_node_dict.values():
    if node not in full_node_dict.values():
        # add node to dict
        new_node_ind = len(full_node_dict)
        full_node_dict.update({new_node_ind:node})
        # add node to full graph
        full_g.add_vertices(1)
        for attr in node.keys():
            full_g.vs[new_node_ind][attr] = node[attr]

#   DEAL WITH EDGES

# get snapshot id
fn = fp.split("/")[-1]
fnsplit = fn.split("_")
ind = fnsplit.index(splitkey)+1
snap = fnsplit[ind]
# add label to edges
g.es['snap'] = snap
start = full_g.ecount()
edge_attr_names = g.es.attributes()
for edge_ind, edge in enumerate(g.es()):
    (local_source, local_target) = edge.tuple
    full_source = full_node_dict.values().index(local_node_dict[local_source])
    full_target = full_node_dict.values().index(local_node_dict[local_target])
    full_g.add_edges([(full_source,full_target)])
    for att in edge_attr_names:
        g.es[start+edge_ind][att] = edge[att]
# # write to gml file
# fnsplit = fn.split("_" + snap + "_")
# out_gml_title = root_directory + order + "/" + fnsplit[0] + "_Union_" + fnsplit[1].split(".")[0][:-2] + order + ".gml"
# print out_gml_title
# comment = "Union of the " + fnsplit[0][:-1] + " files."
# g.write_gml(out_gml_title, comment)
