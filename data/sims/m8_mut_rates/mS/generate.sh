#!/bin/bash
# run on turing.cs.princeton.edu

if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <simulate_executable>" >&2
    exit 1
fi

for mut in {1,10,100}
do
    if [ ! -e $mut ]
    then
        mkdir $mut
    fi
    for s in {0,10,12,2,3,4,5,7,8,9}
    do
        $1 -kP 2 -k 1 -K 5e4 -D 2e-7 -mut ${mut} -m 7 -s $s -P 0.8 -E 0.001 -C 200 -o $mut/ -c ../../coloring.txt
    done
done

