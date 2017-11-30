#!/bin/bash
for m in {m5,m8}
do
    if [ ! -e input_$m ]
    then
        mkdir input_$m
    fi
    for p in {mS,S,M,R}
    do
        for f in ../../../data/sims/$m/$p/reads_seed*.tsv
        do
            s=$(basename $f .tsv | sed -e s/reads_seed//g)
            echo "Generating TreeOmics input for seed $s pattern $p... (m = $m)"
            python convert_reads_to_treeomics.py ../../../data/sims/$m/$p/reads_seed$s.tsv loci_shuffled.txt > input_$m/${p}_seed${s}_variant.txt 2> input_$m/${p}_seed${s}_coverage.txt
        done
    done
done
