#!/bin/bash

echo "m	p	s	method	ARI	RI	recall	precision"
for x in {ancestree,phylowgs,clomial,sciclone,pyclone,pyclone_binomial,pyclone_org,pyclone_org_binomial}
do
    for m in {m5,m8}
    do
        for p in {mS,S,M,R}
        do
            for f in ../../../data/sims/$m/$p/clustering_observed_seed*.txt
            do
                s=$(basename $f .txt | sed -e s/clustering_observed_seed//g)
                ff=$x/$m/$p/$s/clustering.txt
                if [ "$x" == "pyclone_binomial" ]
                then
                    ff=pyclone/$m/$p/$s/clustering_binomial.txt
                fi
                if [ "$x" == "pyclone_org_binomial" ]
                then
                    ff=pyclone_org/$m/$p/$s/clustering_binomial.txt
                fi
                #echo $f $ff
                if [ -e $ff ] && [ -s $ff ]
                then
                    echo -n "$m	$p	$s	$x	"
                    python compare.py $f $ff
                fi
            done
        done
    done
done
