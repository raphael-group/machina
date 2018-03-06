#!/usr/bin/python
import sys

if len(sys.argv) != 2:
    sys.stderr.write("Usage: %s <filename>\n" % sys.argv[0])
    sys.exit(1)

filename = sys.argv[1]

C = set()
with open(filename) as f:
    for line in f:
        s = line.rstrip("\n").split("\t")
        src = s[0].split(';')
        tgt = s[1].split(';')
        
        okSrc = True
        for ss in src:
            okSrc = okSrc and ss.isdigit()

        if okSrc:
            C.add(s[0])
              
        okTgt = True
        for ss in tgt:
            okTgt = okTgt and ss.isdigit()

        if okTgt:
            C.add(s[1])

for c in C:
    print c
