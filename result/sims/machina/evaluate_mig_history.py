#!/usr/bin/python
import sys
import re

def get_mutations(edge_list, u):
    # find parent
    s = re.split("_|;|^", u)
    mutations = set()
    for i in s:
        if i.isdigit():
            mutations.add(int(i))
    #print mutations

    for edge in edge_list:
        uu = edge[0]
        vv = edge[1]
        if vv == u:
            return mutations | get_mutations(edge_list, uu)

    return mutations

def parse_clone_tree(filename_T, filename_l):
    edges = []
    with open(filename_T) as f:
        for line in f:
            s = line.rstrip("\n").split(" ")
            edges += [(s[0], s[1])]

    labeling = {}
    with open(filename_l) as f:
        for line in f:
            s = line.rstrip("\n").split(" ")
            labeling[s[0]] = s[1]

    # find migration edges
    migration_edges = []
    for (u, v) in edges:
        if labeling[u] != labeling[v]:
            migration_edges += [(u,v)]

    return edges, migration_edges

def identify_seeding_clones(edge_list, migration_edge_list):
    res = set()
    for (u,v) in migration_edge_list:
        muts_u = get_mutations(edge_list, u)
        muts_v = get_mutations(edge_list, v)
        res.add(frozenset(muts_u))

    return res

def parse_migration_graph(filename_G):
    edges = []
    with open(filename_G) as f:
        for line in f:
            s = line.rstrip("\n").split(" ")
            edges += [(s[0], s[1])]

    return edges

def multi_graph_to_set(edge_list):
    count = {}
    res = set()
    for edge in edge_list:
        if edge not in count:
            count[edge] = 1
        else:
            count[edge] += 1
        res.add((edge[0], edge[1], count[edge]))
    return res

if __name__ == "__main__":
    if len(sys.argv) != 7:
        sys.stderr.write("Usage: %s <SIMULATED_CLONE_TREE> <SIMULATED_VERTEX_LABELING> <SIMULATED_MIGRATION_GRAPH>"
                         " <INFERRED_CLONE_TREE> <INFERRED_VERTEX_LABELING> <INFERRED_MIGRATION_GRAPH>\n" % sys.argv[0])
        sys.exit(1)

    edges_simulated, mig_edges_simulated = parse_clone_tree(sys.argv[1], sys.argv[2])
    seeding_clones_simulated = identify_seeding_clones(edges_simulated, mig_edges_simulated)

    edges_inferred, mig_edges_inferred = parse_clone_tree(sys.argv[4], sys.argv[5])
    seeding_clones_inferred = identify_seeding_clones(edges_inferred, mig_edges_inferred)

    recall = float(len(seeding_clones_inferred & seeding_clones_simulated)) / float(len(seeding_clones_simulated))
    precision = float(len(seeding_clones_inferred & seeding_clones_simulated)) / float(len(seeding_clones_inferred))
    if recall == 0 or precision == 0:
        F = 0
    else:
        F = 2.0 / ((1.0 / recall) + (1.0 / precision))

    edge_set_G_simulated = set(parse_migration_graph(sys.argv[3]))
    edge_set_G_inferred = set(parse_migration_graph(sys.argv[6]))

    edge_multiset_G_simulated = multi_graph_to_set(parse_migration_graph(sys.argv[3]))
    edge_multiset_G_inferred = multi_graph_to_set(parse_migration_graph(sys.argv[6]))

    recall_G = float(len(edge_set_G_inferred & edge_set_G_simulated)) / float(len(edge_set_G_simulated))
    precision_G = float(len(edge_set_G_inferred & edge_set_G_simulated)) / float(len(edge_set_G_inferred))

    recall_G2 = float(len(edge_multiset_G_inferred & edge_multiset_G_simulated)) / float(len(edge_multiset_G_simulated))
    precision_G2 = float(len(edge_multiset_G_inferred & edge_multiset_G_simulated)) / float(len(edge_multiset_G_inferred))

    if recall_G != 0 and precision_G != 0:
        F_G = 2.0 / ((1.0 / recall_G) + (1.0 / precision_G))
    else:
        F_G = 0
    if recall_G2 != 0 and precision_G2 != 0:
        F_G2 = 2.0 / ((1.0 / recall_G2) + (1.0 / precision_G2))
    else:
        F_G2 = 0

    print ",".join(map(str, [recall, precision, F, recall_G, precision_G, F_G, recall_G2, precision_G2, F_G2]))
