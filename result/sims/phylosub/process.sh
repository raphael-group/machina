#!/bin/bash
if [ ! $# -eq 2 ]
then
    echo "Usage: $0 <clonetreevisualization_executable> <rf_executable>" >&2
    exit 1
fi

for m in {m5,m8}
do
    echo "pattern,seed,method,RF" > results_${m}.txt
    for p in {mS,S,M,R}
    do
        for f in ../../../data/sims/${m}/$p/reads_seed*.tsv
        do
	    for filter in {0,1,2,5,10}
	    do
		s=$(basename $f .tsv | sed -e s/reads_seed//g)

		echo Solving seed $s, pattern $p, anatomical sites $m...
		$1 -c ../../../data/sims/coloring.txt ${m}/${p}/reads_seed${s}.phylowgs/edgelist_${filter}.txt ${m}/${p}/reads_seed${s}.phylowgs/leaflabels_${filter}.txt > ${m}/${p}/reads_seed${s}.phylowgs/${filter}_tree.dot
		$2 ../../../data/sims/${m}/$p/T_seed$s.tree ../../../data/sims/${m}/$p/T_seed$s.labeling ${m}/${p}/reads_seed${s}.phylowgs/edgelist_${filter}.txt ${m}/${p}/reads_seed${s}.phylowgs/leaflabels_${filter}.txt > ${m}/${p}/reads_seed${s}.phylowgs/RF_${filter}.txt
		echo -n $p,$s,PhyloWGS $filter, >> results_${m}.txt
		tail -n 1 ${m}/${p}/reads_seed${s}.phylowgs/RF_${filter}.txt | cut -d' ' -f3  >> results_${m}.txt
	    done
        done
    done
done
