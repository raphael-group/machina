#!/bin/bash
if [ ! $# -eq 2 ]
then
    echo "Usage: $0 <cluster> <m>" >&2
    exit 1
fi

m=$2

if [ ! -e input_${m} ]
then
    mkdir input_${m}
fi


for p in {mS,S,M,R}
do
    for f in ../../../data/sims/${m}/$p/reads_seed*.tsv
    do
	s=$(basename $f .tsv | sed -e s/reads_seed//g)
	echo Clustering seed $s pattern $p...
	$1 -a 0.001 -b 0.05 $f > input_${m}/cluster_${p}_seed${s}.tsv 2> input_${m}/cluster_${p}_seed${s}.txt
    done
done
