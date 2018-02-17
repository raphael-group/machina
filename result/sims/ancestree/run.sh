#!/bin/bash
if [ ! $# -eq 5 ]
then
    echo "Usage: $0 <cluster_executable> <ancestree_executable> <machina_input_executable> <clonetreevisualization_executable> <rf_executable>" >&2
    exit 1
fi

for m in {m5,m8}
do
    if [ ! -e input_${m} ]
    then
	mkdir input_${m}
    fi

    if [ ! -e result_${m} ]
    then
	mkdir result_${m}
    fi

    echo "pattern,seed,method,RF" > results_${m}.txt
    for p in {mS,S,M,R}
    do
        for f in ../../../data/sims/${m}/$p/reads_seed*.tsv
        do
    	    s=$(basename $f .tsv | sed -e s/reads_seed//g)

	    echo Solving seed $s, pattern $p, anatomical sites $m...
	    #$1 -a 0.001 -b 0.05 -A $f > input_${m}/cluster_${p}_seed${s}.tsv 2> input_${m}/cluster_${p}_seed${s}.txt
	    #$2 input_${m}/cluster_${p}_seed${s}.tsv 2> result_${m}/${p}_seed${s}.log | $3 0 - 0.05 0.8 > result_${m}/${p}_seed${s}.tree 2> result_${m}/${p}_seed${s}.labeling
	    #sed -i '' -e 's/^\(.*_\(.*\)_.*\)$/\1 \2/g' result_${m}/${p}_seed${s}.labeling
	    $4 -c ../../../data/sims/coloring.txt result_${m}/${p}_seed${s}.tree result_${m}/${p}_seed${s}.labeling > result_${m}/${p}_seed${s}.dot
	    $5 ../../../data/sims/${m}/$p/T_seed$s.tree ../../../data/sims/${m}/$p/T_seed$s.labeling result_${m}/${p}_seed${s}.tree result_${m}/${p}_seed${s}.labeling > result_${m}/${p}_seed${s}.RF.txt
	    echo -n $p,$s,AncesTree, >> results_${m}.txt
	    tail -n 1 result_${m}/${p}_seed${s}.RF.txt | cut -d' ' -f3  >> results_${m}.txt
        done
    done
done
