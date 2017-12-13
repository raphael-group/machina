#!/bin/bash
for m in {m5,m8}
do
    for x in {_pyclone,_sciclone,_phylowgs,_clomial,""}
    do
        echo $m $x
        ./process_pmh_cti.sh output_${m}${x}/ ../../../build/rf $m > results_MACHINA${x}_${m}.txt
    done
done
