#!/bin/bash
for m in {m5,m8}
do
    for p in {mS,S,M,R}
    do
        for f in ../../../data/sims/$m/$p/reads_seed*.tsv
        do
            s=$(basename $f .tsv | sed -e s/reads_seed//g)
            echo Running pyclone input for $m, seed $s, pattern $p...
	    python ./parse_pyclone_output.py pyclone/$m/$p/$s/clustering.tsv > pyclone/$m/$p/$s/clustering.txt
	    python ./parse_pyclone_output.py pyclone/$m/$p/$s/clustering_binomial.tsv > pyclone/$m/$p/$s/clustering_binomial.txt
        done
    done
done
