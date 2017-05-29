#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <sankoff_executable>" >&2
    exit 1
fi

if [ ! -d sankoff ]
then
    mkdir sankoff
fi

for i in {1,2,3,4,7,9,10}
do
    echo "Solving patient${i}..."
    if [ ! -d sankoff/patient${i} ]
    then
        mkdir sankoff/patient${i}
    fi
    $1 ../../data/mcpherson_2016/patient${i}.tree ../../data/mcpherson_2016/patient${i}.labeling -p LOv,ROv,RUt -c ../../data/mcpherson_2016/coloring.txt -o sankoff/patient${i}/ 2> sankoff/patient${i}/result.txt
done
