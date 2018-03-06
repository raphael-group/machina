#!/bin/bash
if [ ! $# -eq 2 ]
then
    echo "Usage: $0 <ms_executable> <pmh_tr_executable>" >&2
    exit 1
fi

#for p in {A10,A22,A29,A31,A32}
for p in {A10,A31,A32}
do
    echo Running $p...

    if [ ! -e $p ]
    then
	mkdir $p
    fi

    #f=../../data/gundem_2015/$p.tsv
    #$1 -p prostate $f 2> $p/result_ms.txt

    #python plot_frequency_intervals.py $f 9 colors.txt $p/`basename $f .tsv`.pdf

    t=../../data/gundem_2015/reported_clonetrees/$p.tree
    l=../../data/gundem_2015/reported_clonetrees/$p.labeling

    $2 -p prostate -c ../../data/gundem_2015/coloring.txt -o $p $t $l > $p/result_pmh-tr.txt
done
