#!/bin/bash
if [ ! $# -eq 2 ]
then
    echo "Usage: $0 <generatemutationtrees> <m>" >&2
    exit 1
fi

m=$2

if [ ! -e mut_trees_${m} ]
then
    mkdir mut_trees_${m}
fi

for p in {mS,S,M,R}
do
    for f in ../../../data/sims/${m}/$p/reads_seed*.tsv
    do
	s=$(basename $f .tsv | sed -e s/reads_seed//g)
	echo Generating mutation trees: seed $s pattern $p...
	$1 input_${m}/cluster_${p}_seed${s}.tsv > mut_trees_${m}/mut_trees_${p}_seed${s}.txt
    done
done
