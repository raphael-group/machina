#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <pmh_pr_executable>" >&2
    exit 1
fi

if [ ! -d pmh_pr ]
then
    mkdir pmh_pr
fi

for i in {1,2,3,4,7,9,10}
do
    echo "Solving patient${i}..."
    if [ ! -d pmh_pr/patient${i} ]
    then
	mkdir pmh_pr/patient${i}
    fi
    $1 ../../data/mcpherson_2016/patient${i}.tree ../../data/mcpherson_2016/patient${i}.labeling -t 4 -l 3600 -p LOv,ROv,RUt -c ../../data/mcpherson_2016/coloring.txt -o pmh_pr/patient${i}/ 2> pmh_pr/patient${i}/result.txt
done
#wait
