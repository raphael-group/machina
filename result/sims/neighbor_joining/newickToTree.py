#!/bin/bash
from Bio import Phylo
import sys

def tabulate_names(tree):
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if clade.name:
            pass
            #clade.name = '%d_%s' % (idx, clade.name)
        else:
            clade.name = str(idx)
        names[clade.name] = clade
    return names

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: %s <NEWICK>\n" % sys.argv[0])
        sys.exit(1)

    tree = Phylo.read(sys.argv[1], "newick")
    tree.rooted = True
    table = tabulate_names(tree)
    n = Phylo.to_networkx(tree)

    print "GL", 0
    for edge in n.edges():
        print edge[0], edge[1]

    with open(sys.argv[1] + ".ascii", "w") as f:
        Phylo.draw_ascii(tree, file=f)

    for t in tree.get_terminals():
        s = str(t).split("_")
        sys.stderr.write("%s %s\n" % (t, s[0]))
