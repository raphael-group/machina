#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <pmh_sankoff_executable>" >&2
    exit 1
fi

if [ ! -d pmh_sankoff ]
then
    mkdir pmh_sankoff
fi

for i in {1,2,3,4,7,9,10}
do
    echo "Solving patient${i}..."
    if [ ! -d pmh_sankoff/patient${i} ]
    then
        mkdir pmh_sankoff/patient${i}
    fi
    $1 ../../data/mcpherson_2016/patient${i}.tree ../../data/mcpherson_2016/patient${i}.labeling -p LOv,ROv,RUt -c ../../data/mcpherson_2016/coloring.txt -o pmh_sankoff/patient${i}/ 2> pmh_sankoff/patient${i}/result.txt
done
