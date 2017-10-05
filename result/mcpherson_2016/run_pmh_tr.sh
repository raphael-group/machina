#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <pmh_tr_executable>" >&2
    exit 1
fi

if [ ! -d pmh_tr ]
then
    mkdir pmh_tr
fi

for i in {1,2,3,4,7,9,10}
do
    echo "Solving patient${i}..."
    if [ ! -d pmh_tr/patient${i} ]
    then
	mkdir pmh_tr/patient${i}
    fi
    $1 ../../data/mcpherson_2016/patient${i}.tree ../../data/mcpherson_2016/patient${i}.labeling -t 4 -p LOv -c ../../data/mcpherson_2016/coloring.txt -o pmh_tr/patient${i}/ > pmh_tr/patient${i}/result_LOv.txt
    $1 ../../data/mcpherson_2016/patient${i}.tree ../../data/mcpherson_2016/patient${i}.labeling -t 4 -p ROv -c ../../data/mcpherson_2016/coloring.txt -o pmh_tr/patient${i}/ > pmh_tr/patient${i}/result_ROv.txt
done
