#!/bin/bash
if [ ! $# -eq 2 ]
then
    echo "Usage: $0 <pmh_ti_executable> <generatemutationtrees_executable>" >&2
    exit 1
fi

if [ ! -e A1 ]
then
    mkdir A1
fi

$2 ../../data/hoadley_2016/A1/A1_MACHINA_0.95.tsv > A1/mutation_trees.txt

$1 -c ../../data/hoadley_2016/coloring.txt -barT A1/mutation_trees.txt -F ../../data/hoadley_2016/A1/A1_MACHINA_0.95.tsv -p breast -t 4 -o A1/ > A1/results.txt
