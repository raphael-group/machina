#!/bin/bash
if [ ! $# -eq 3 ]
then
    echo "Usage: $0 <pmh_ti_executable> <generatemigrationtrees_executable> <generatemutationtrees_executable>" >&2
    exit 1
fi

if [ ! -e A7 ]
then
    mkdir A7
fi

$2 breast brain kidney liver lung rib > A7/migration_trees.txt

$3 ../../data/hoadley_2016/A7/F.tsv > A7/mutation_trees.txt

$1 -c ../../data/hoadley_2016/coloring.txt -barT A7/mutation_trees.txt -G A7/migration_trees.txt -F ../../data/hoadley_2016/A7/A7_MACHINA_0.95.tsv -p breast -t 4 -m 1 -o A7/ > A7/results.txt
