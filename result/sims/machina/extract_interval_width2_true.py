#!/usr/bin/python
import sys
import os

if len(sys.argv) != 3:
    sys.stderr.write("Usage: %s <FREQ_INPUT_FILE> <METHOD>\n" % sys.argv[0])
    sys.exit(1)

#R/clustering_observed_seed247.txt
filename = sys.argv[1]
instance = os.path.basename(filename).lstrip("reads_clustering_observed_").rstrip(".tsv")
pattern = filename.split("/")[-2]
instance = instance.split("_")[0].lstrip("seed")
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
        print "\t".join(map(str, [sys.argv[2],instance, sample, character, width]))
