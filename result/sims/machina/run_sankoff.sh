#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <pmh_sankoff_executable>" >&2
    exit 1
fi

if [ ! -e sankoff ]
then
    mkdir sankoff
fi

for m in {m5,m8}
do
    for p in {mS,S,M,R}
    do
        for f in ../../../data/sims/$m/$p/reads_seed*.tsv
        do
            s=$(basename $f .tsv | sed -e s/reads_seed//g)
            echo Solving seed $s, pattern $p, anatomical sites $m...
            tree=../../../data/sims/$m/$p/T_seed${s}.tree
            labeling=../../../data/sims/$m/$p/T_seed${s}.labeling
            $1 -p P $tree $labeling 2> sankoff/${m}_${p}_${s}.txt
        done
    done
done
