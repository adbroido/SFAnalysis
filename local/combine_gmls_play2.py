import os
import glob
import igraph
import numpy as np

root_directory = "/Users/annabroido/CU/testgmls/Social/Political/"
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

full_g = igraph.read(paths[0])
partial_g = igraph.read(paths[1])

full_nodes = [node.attributes() for node in full_g.vs]
partial_nodes = [node.attributes() for node in partial_g.vs]

import itertools
new_nodes = list(itertools.ifilterfalse(lambda x: x in full_nodes, partial_nodes))
