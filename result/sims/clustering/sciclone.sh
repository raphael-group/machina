#!/bin/bash
if [ ! -e sciclone ]
then
    mkdir sciclone
fi

for m in {m5,m8}
do
    if [ ! -e sciclone/$m ]
    then
        mkdir sciclone/$m
    fi
    for p in {mS,S,M,R}
    do
        if [ ! -e sciclone/$m/$p ]
        then
            mkdir sciclone/$m/$p
        fi
        for f in ../../../data/sims/${m}/$p/reads_seed*.tsv
        do
            s=$(basename $f .tsv | sed -e s/reads_seed//g)
            if [ ! -e sciclone/$m/$p/$s ]
            then
                mkdir sciclone/$m/$p/$s
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
            echo Running sciclone input for $m, seed $s, pattern $p...
            #python ./convert_reads_to_sciclone.py ../../../data/sims/$m/$p/reads_seed${s}.tsv sciclone/$m/$p/$s/
            #Rscript sciclone.R sciclone/${m}/${p}/$s 10
	    python ./parse_sciclone_output.py sciclone/$m/$p/$s/clustering.tsv > sciclone/$m/$p/$s/clustering.txt
        done
    done
done
