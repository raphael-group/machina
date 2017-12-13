#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <simulate_executable>" >&2
    exit 1
fi

#INITIAL run, which generated 10 random seeds resulting in pS
$1 -kP 3 -k 3 -D 2e-7 -m 4 -p 1 -s 0 -N 10 -C 10000 -o . -c ../../coloring.txt

#for s in {10,15,26,27,42,52,54,66,8,9}
#do
#    $1 -kP 5 -k 5 -D 2e-7 -m 4 -p 1 -s $s -P 1 -C 10000 -o . -c ../../coloring.txt
#done
