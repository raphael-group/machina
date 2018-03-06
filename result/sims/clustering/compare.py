#!/usr/bin/python
import clustering as clust
import sys

def parse_solution(filename):
    C = []
    with open(filename) as f:
        for line in f:
            s = line.rstrip("\n").split(';')
            cc = []
            for ss in s:
                if ss.isdigit():
                    cc.append(int(ss))
            C += [cc]

    return C

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: %s <true_solution> <inferred_solution>\n" % sys.argv[0])
        sys.exit(1)

    true_sol_filename = sys.argv[1]
    inf_sol_filename = sys.argv[2]

    C_true = parse_solution(true_sol_filename)
    C_inf = parse_solution(inf_sol_filename)

    ari = clust.adjusted_rand_index(C_true, C_inf)
    ri = clust.rand_index(C_true, C_inf)
    recall, precision = clust.recall_and_precision(C_inf, C_true)

    print "\t".join([str(ari), str(ri), str(recall), str(precision)])
