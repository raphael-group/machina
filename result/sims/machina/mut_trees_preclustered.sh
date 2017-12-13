#!/bin/bash
if [ ! $# -eq 4 ]
then
    echo "Usage: $0 <cluster> <generatemutationtrees> <m> <clustering_dir>" >&2
    exit 1
fi

m=$3
x=$4
xx=$(basename $x)

if [ ! -e mut_trees_${m}_${xx} ]
then
    mkdir mut_trees_${m}_${xx}
fi

for p in {mS,S,M,R}
do
    for f in ../../../data/sims/${m}/$p/reads_seed*.tsv
    do
	s=$(basename $f .tsv | sed -e s/reads_seed//g)
        if [ $xx == "pyclone_org" ]
        then
            pre="${x}/${m}/${p}/${s}/clustering_binomial.txt"
        else
            pre="${x}/${m}/${p}/${s}/clustering.txt"
        fi
	echo Generating mutation trees: seed $s pattern $p...
        $1 -a 0.001 -b 0.05 -C ${pre} $f > mut_trees_${m}_${xx}/cluster_${p}_seed${s}.tsv 2> mut_trees_${m}_${xx}/cluster_${p}_seed${s}.txt
	$2 mut_trees_${m}_${xx}/cluster_${p}_seed${s}.tsv > mut_trees_${m}_${xx}/mut_trees_${p}_seed${s}.txt
    done
done
