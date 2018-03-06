#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <simulate_executable>" >&2
    exit 1
fi

#INITIAL run, which generated 10 random seeds resulting in mS
$1 -kP 3 -k 3 -D 2e-7 -m 4 -p 0 -s 0 -N 10 -C 10000 -o . -c ../../coloring.txt

#subsequent runs
#for s in {0,10,12,2,3,4,5,7,8,9}
#do
#    $1 -kP 5 -k 5 -D 2e-7 -m 4 -p 0 -s $s -P 1 -C 10000 -o . -c ../../coloring.txt
#done
