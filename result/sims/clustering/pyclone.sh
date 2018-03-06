#!/bin/bash
if [ ! -e pyclone ]
then
    mkdir pyclone
fi

for m in {m5,m8}
do
    if [ ! -e pyclone/$m ]
    then
        mkdir pyclone/$m
    fi
    for p in {mS,S,M,R}
    do
        if [ ! -e pyclone/$m/$p ]
        then
            mkdir pyclone/$m/$p
        fi
        for f in ../../../data/sims/${m}/$p/reads_seed*.tsv
        do
            s=$(basename $f .tsv | sed -e s/reads_seed//g)
            if [ ! -e pyclone/$m/$p/$s ]
            then
                mkdir pyclone/$m/$p/$s
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
            echo Generating PyClone input for $m, seed $s, pattern $p...
            python ./convert_reads_to_pyclone.py ../../../data/sims/$m/$p/reads_seed${s}.tsv pyclone/$m/$p/$s/
    
            d=pyclone/$m/$p/$s
            cp config.yaml $d/
            cp config_binomial.yaml $d/config_binomial.yaml
    
            for f in pyclone/$m/$p/$s/*.tsv
            do
                if [ `basename $f` == "clustering.tsv" ]; then continue; fi
                if [ `basename $f` == "clustering_binomial.tsv" ]; then continue; fi
                ff=`dirname $f`/`basename $f .tsv`.yaml
                #if [ ! -e $ff ]
                #then
                PyClone build_mutations_file --in_file $f --out_file $ff &
                #fi
                echo >> $d/config.yaml
                echo "  `basename $f .tsv`:" >> $d/config.yaml
                echo "    mutations_file: `basename $ff`" >> $d/config.yaml
                echo "    tumour_content:" >> $d/config.yaml
                echo "      value: 1.0" >> $d/config.yaml
                echo "    error_rate: 0.001" >> $d/config.yaml
    
                echo >> $d/config_binomial.yaml
                echo "  `basename $f .tsv`:" >> $d/config_binomial.yaml
                echo "    mutations_file: `basename $ff`" >> $d/config_binomial.yaml
                echo "    tumour_content:" >> $d/config_binomial.yaml
                echo "      value: 1.0" >> $d/config_binomial.yaml
                echo "    error_rate: 0.001" >> $d/config_binomial.yaml
    
            done
            wait
        done
    done
done
