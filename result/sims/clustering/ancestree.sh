#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <cluster_executable>" >&2
    exit 1
fi

if [ ! -e ancestree ]
then
    mkdir ancestree
fi

for m in {m5,m8}
do
    if [ ! -e ancestree/${m} ]
    then
	mkdir ancestree/${m}
    fi

    for p in {mS,S,M,R}
    do
        if [ ! -e ancestree/$m/$p ]
        then
            mkdir ancestree/$m/$p
        fi

        for f in ../../../data/sims/${m}/$p/reads_seed*.tsv
        do
    	    s=$(basename $f .tsv | sed -e s/reads_seed//g)
            if [ ! -e ancestree/$m/$p/$s ]
            then
                mkdir ancestree/$m/$p/$s
            fi

	    echo Solving seed $s, pattern $p, anatomical sites $m...
	    $1 -a 0.001 -b 0.05 $f 2> ancestree/${m}/${p}/${s}/clustering.txt > /dev/null
        done
    done
done
