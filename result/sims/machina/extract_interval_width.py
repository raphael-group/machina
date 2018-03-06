#!/usr/bin/python
import sys
import os

if len(sys.argv) != 2:
    sys.stderr.write("Usage: %s <FREQ_INPUT_FILE>\n" % sys.argv[0])
    sys.exit(1)

filename = sys.argv[1]
instance = os.path.basename(filename).lstrip("cluster_").rstrip(".tsv")
pattern = instance.split("_")[0]
instance = instance.split("_")[1].lstrip("seed")
#print "\t".join(["method", "instance", "sample", "character", "width"])
L = []
with open(filename) as f:
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    for line in f:
        s = line.rstrip("\n").split("\t")
        sample = s[0]
        character = s[4]
        lb = float(s[6])
        ub = float(s[7])
        width = ub - lb
        L.append(width)

import numpy as np
print np.mean(L),
