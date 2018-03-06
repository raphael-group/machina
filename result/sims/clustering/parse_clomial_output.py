#!/usr/bin/python
import sys

def get_nrClusters(dirname):
    min_bic = float('inf')
    nrClusters = -1
    with open(dirname + "/" + "bics.tsv") as f:
        f.readline()
        for line in f:
            s = line.replace('"', '').rstrip("\n").split("\t")
            if s[1] == "NA":
                continue
            bic = float(s[1])
            if bic < min_bic:
                min_bic = bic
                nrClusters = int(s[0])

    return nrClusters

if len(sys.argv) != 2:
    sys.stderr.write("Usage: %s <dirname>\n" % sys.argv[0])
    sys.exit(1)

dirname = sys.argv[1]
nrClusters = get_nrClusters(dirname)

with open(dirname + "/" + str(nrClusters) + ".tsv") as f:
    f.readline()
    C = {}
    for line in f:
        s = line.replace('"', '').rstrip("\n").split("\t")
        pos = s[0].split(":")[1]
        genotype = "".join(s[1:])
        if genotype not in C:
            C[genotype] = []
        C[genotype].append(pos)

    for genotype in C:
        print ";".join(C[genotype])
