#!/bin/bash
if [ ! -e phylowgs ]
then
    mkdir phylowgs
fi

for m in {m5,m8}
do
    if [ ! -e phylowgs/$m ]
    then
        mkdir phylowgs/$m
    fi
    for p in {mS,S,M,R}
    do
        if [ ! -e phylowgs/$m/$p ]
        then
            mkdir phylowgs/$m/$p
        fi
        for f in ../../../data/sims/${m}/$p/reads_seed*.tsv
        do
            s=$(basename $f .tsv | sed -e s/reads_seed//g)
            if [ ! -e phylowgs/$m/$p/$s ]
            then
                mkdir phylowgs/$m/$p/$s
            fi
        done
    done
done

for m in {m5,m8}
do
    for p in {mS,S,M,R}
    do
        for f in ../../../data/sims/$m/$p/reads_seed*.tsv
        do
            s=$(basename $f .tsv | sed -e s/reads_seed//g)
            echo Running phylowgs input for $m, seed $s, pattern $p...
            #python ./convert_reads_to_phylowgs.py ../../../data/sims/$m/$p/reads_seed${s}.tsv phylowgs/$m/$p/$s/
            #Rscript phylowgs.R phylowgs/${m}/${p}/$s 10
	    python ./parse_phylowgs_output.py ../phylosub/$m/$p/reads_seed${s}.phylowgs/edgelist.txt > phylowgs/$m/$p/$s/clustering.txt
        done
    done
done
