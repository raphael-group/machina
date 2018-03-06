#!/bin/bash

echo "method	pattern	instance	sample	character	width"
for f in ../../../data/sims/m5/*/reads_clustering*.tsv
do
    python extract_interval_width2.py $f TRUE
done

for method in {clomial,phylowgs,pyclone,sciclone,ancestree}
do
    for f in mut_trees_m5_${method}/cluster*.tsv
    do
        python extract_interval_width2.py $f $method
    done
done
