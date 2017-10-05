#!/bin/bash
if [ ! $# -eq 3 ]
then
    echo "Usage: $0 <cluster_executable> <generatemutationtrees_executable> <pmh_ti_executable>" >&2
    exit 1
fi

p=A22_noCNA
if [ ! -e $p ]
then
    mkdir $p
fi

f=../../data/gundem_2015/A22_noCNA/reads_A22.tsv

$1 -r $f > $p/frequencies_A22.tsv 2> $p/clusters_A22.tsv
$2 $p/frequencies_A22.tsv > $p/mutationtrees_A22.txt
$3 -p prostate -c ../../data/gundem_2015/coloring.txt -F $p/frequencies_A22.tsv -barT $p/mutationtrees_A22.txt -o $p
