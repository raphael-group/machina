#!/bin/bash
# run on turing.cs.princeton.edu

if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <downsample_executable>" >&2
    exit 1
fi

for p in {mS,S,M,R}
do
    for f in ../m8_mut_rates/$p/reads_seed*_mut10.tsv
    do
        for s in {0.005,0.01,0.05,0.1,0.5,1.0}
        do
            $1 -s $RANDOM -snvFrac $s $f > $p/`basename $f .tsv`_SNV${s}.tsv
        done
    done
done
