#!/bin/bash
if [ ! $# -eq 3 ]
then
    echo "Usage: $0 <cluster_executable> <generatemutationtrees_executable> <pmh_ti_executable>" >&2
    exit 1
fi

if [ ! -e input ]
then
    mkdir input
fi

if [ ! -e mut_trees ]
then
    mkdir mut_trees
fi

if [ ! -e output ]
then
    mkdir output
fi

for f in ../../data/sanborn_2015/reads_*.tsv
do
    s=$(basename $f .tsv | sed -e s/reads_//g)
    echo Clustering seed `basename $f`...
    if [ $s == "B" ]
    then
	$1 -r -a 0.001 -b 0.001 $f > input/cluster_${s}.tsv 2> input/cluster_${s}.txt
    else
	$1 -r -a 0.001 $f > input/cluster_${s}.tsv 2> input/cluster_${s}.txt
    fi
done

for f in ../../data/sanborn_2015/reads_*.tsv
do
    s=$(basename $f .tsv | sed -e s/reads_//g)
    echo Generating mutation trees: $s...
    $2 input/cluster_${s}.tsv > mut_trees/mut_trees_${s}.txt
done

for p in {A,B,C,D,E,F,G}
do
    echo Running patient $p...
    if [ ! -e output/$p ]
    then
	mkdir output/$p
    fi
    $3 -c ../../data/sanborn_2015/coloring.txt -barT mut_trees/mut_trees_${p}.txt -F input/cluster_${p}.tsv -o output/$p -p primary > output/result_$p.txt
done
