#!/bin/bash
if [ ! $# -eq 2 ]
then
    echo "Usage: $0 <pmh_cti_executable> <generatemigrationtrees_executable>" >&2
    exit 1
fi

if [ ! -e A7 ]
then
    mkdir A7
fi

if [ ! -e A7/sample_trees ]
then
    mkdir A7/sample_trees
fi
cd A7/sample_trees
../../$2 breast brain kidney liver lung rib
cd ../..

for T in A7/sample_trees/*.txt
do
    for i in {0,1}
    do
	dir=A7/T`basename ${T} .txt`_sol${i}
	if [ ! -e $dir ]
	then
	    mkdir $dir
	fi
	echo "Solving $dir..."
	$1 -c ../../data/hoadley_2016/coloring.txt ../../data/hoadley_2016/A7/sol${i}.tree ../../data/hoadley_2016/A7/F.tsv -p breast -t 4 -m 1 -o $dir -l 30 -s $T &> $dir/result.txt
    done
done
