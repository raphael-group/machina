#!/usr/bin/python
import sys

with open(sys.argv[1]) as f:
    print f.readline(),
    print f.readline(),
    for line in f:
        s = line.split("\t")
        s[3] = "(" + s[3] + ",(1,1,1))"
        print "\t".join(s[0:5] + [s[6]])
