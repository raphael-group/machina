#!/usr/bin/python
import sys

if len(sys.argv) != 2:
    sys.stderr.write("Usage: %s <clustering.tsv>\n" % sys.argv[0])
    sys.exit(1)

filename = sys.argv[1]

C = {}
with open(filename) as f:
    f.readline()
    for line in f:
        s = line.rstrip("\n").split("\t")
        snv = str(int(s[0]) / 1000)
        cluster = s[2]
        if cluster not in C:
           C[cluster] = set()
        C[cluster].add(snv)

    for cluster in C:
        print ";".join(list(C[cluster]))

