#!/bin/bash
if [ ! -e clomial ]
then
    mkdir clomial
fi

for m in {m5,m8}
do
    if [ ! -e clomial/$m ]
    then
        mkdir clomial/$m
    fi
    for p in {mS,S,M,R}
    do
        if [ ! -e clomial/$m/$p ]
        then
            mkdir clomial/$m/$p
        fi
        for f in ../../../data/sims/${m}/$p/reads_seed*.tsv
        do
            s=$(basename $f .tsv | sed -e s/reads_seed//g)
            if [ ! -e clomial/$m/$p/$s ]
            then
                mkdir clomial/$m/$p/$s
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
            echo Running clomial input for $m, seed $s, pattern $p...
            #python ./convert_reads_to_clomial.py ../../../data/sims/$m/$p/reads_seed${s}.tsv clomial/$m/$p/$s/
            #Rscript clomial.R clomial/${m}/${p}/$s 10
	    python ./parse_clomial_output.py clomial/$m/$p/$s > clomial/$m/$p/$s/clustering.txt
        done
    done
done
