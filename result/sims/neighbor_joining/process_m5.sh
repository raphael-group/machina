#!/bin/bash
if [ ! $# -eq 2 ]
then
    echo "Usage: $0 <clonetreevisualization_executable> <rf_executable>" >&2
    exit 1
fi

if [ ! -e result_m5 ]
then
    mkdir result_m5
fi

echo "pattern,seed,method,RF" > results_m5.txt
for p in {mS,S,M,R}
do
    for f in ../../../data/sims/m5/$p/reads_seed*.tsv
    do
	s=$(basename $f .tsv | sed -e s/reads_seed//g)
	echo Running neighbor joining for seed $s pattern $p...
	if [ ! -e ${p}_seed${s}.txt ]
	then
	    python convert_reads_to_nj.py ../../../data/sims/m5/$p/reads_seed$s.tsv > result_m5/${p}_seed${s}.txt
	fi
	if [ ! -e result_m5/${p}_seed${s}.newick ]
	then
	    Rscript run_nj.R result_m5/${p}_seed${s}.txt > result_m5/${p}_seed${s}.newick
	fi
	if [ ! -e ${p}_seed${s}.tree ]
	then
	    python newickToTree.py result_m5/${p}_seed${s}.newick > result_m5/${p}_seed${s}.tree 2> result_m5/${p}_seed${s}.labeling
	fi
	if [ ! -e ${p}_seed${s}.dot ]
	then
	    $1 -c ../../../data/sims/coloring.txt result_m5/${p}_seed${s}.tree result_m5/${p}_seed${s}.labeling > result_m5/${p}_seed${s}.dot
	fi
	$2 ../../../data/sims/m5/$p/T_seed$s.tree ../../../data/sims/m5/$p/T_seed$s.labeling result_m5/${p}_seed${s}.tree result_m5/${p}_seed${s}.labeling > result_m5/${p}_seed${s}.RF.txt
	echo -n $p,$s,Neighbor joining, >> results_m5.txt
       	tail -n 1 result_m5/${p}_seed${s}.RF.txt | cut -d' ' -f3  >> results_m5.txt
    done
done
