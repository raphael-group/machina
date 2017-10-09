#!/bin/bash
if [ ! $# -eq 2 ]
then
    echo "Usage: $0 <clonetreevisualization_executable> <rf_executable>" >&2
    exit 1
fi

echo "pattern,seed,method,RF" > results_m5_sub.txt
for p in {mS,S,M,R}
do
    for f in ../../../data/sims/m5/$p/reads_seed*.tsv
    do
        s=$(basename $f .tsv | sed -e s/reads_seed//g)
        echo Processing TreeOmics output for seed $s pattern $p...
        if [ -e /n/fs/ragr-data/projects/MACHINA/treeomics/src/m5/outputSUB/${p}_seed${s}SUBCLONES/*mlhtree_1 ]
        then
            python transform.py /n/fs/ragr-data/projects/MACHINA/treeomics/src/m5/outputSUB/${p}_seed${s}SUBCLONES/*mlhtree_1 > result_m5_sub/${p}_seed${s}.tree 2> result_m5_sub/${p}_seed${s}.labeling
        elif [ -e /n/fs/ragr-data/projects/MACHINA/treeomics/src/m5/outputSUB/${p}_seed${s}SUBCLONES/*mlhtree_2 ]
        then
            python transform.py /n/fs/ragr-data/projects/MACHINA/treeomics/src/m5/outputSUB/${p}_seed${s}SUBCLONES/*mlhtree_2 > result_m5_sub/${p}_seed${s}.tree 2> result_m5_sub/${p}_seed${s}.labeling
        fi
        $1 -c ../../../data/sims/coloring.txt result_m5_sub/${p}_seed${s}.tree result_m5_sub/${p}_seed${s}.labeling > result_m5_sub/${p}_seed${s}.dot
        $2 ../../../data/sims/m5/$p/T_seed$s.tree ../../../data/sims/m5/$p/T_seed$s.labeling result_m5_sub/${p}_seed${s}.tree result_m5_sub/${p}_seed${s}.labeling > result_m5_sub/${p}_seed${s}.RF.txt
        echo -n $p,$s,Treeomics, >> results_m5_sub.txt
        tail -n 1 result_m5_sub/${p}_seed${s}.RF.txt | cut -d' ' -f3  >> results_m5_sub.txt
    done
done
