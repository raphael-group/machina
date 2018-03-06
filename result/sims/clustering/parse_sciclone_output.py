#!/usr/bin/python
import sys
import pandas as pd

if len(sys.argv) != 2:
    sys.stderr.write("Usage: %s <clustering.tsv>\n" % sys.argv[0])
    sys.exit(1)

filename = sys.argv[1]
df = pd.read_csv(filename, sep="\t")

nrClusters = int(max(df['cluster']))
nrSNVs = max(df['st']) / 1000 + 1

C = [ [] for j in range(nrClusters) ]

for i in range(nrSNVs):
    if len(df[df['st'] == i * 1000]) == 1:
        j = int(df[df['st'] == i * 1000]['cluster'])
        C[j - 1].append(i)

for j in range(nrClusters):
    print ";".join(map(str, C[j]))
