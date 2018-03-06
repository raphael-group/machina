#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <simulate_executable>" >&2
    exit 1
fi

$1 -kP 3 -k 3 -D 2e-7 -m 4 -p 2 -s 0 -N 10 -C 10000 -o . -c ../../coloring.txt

#for s in {0,114,139,191,24,249,250,4,54,98}
#do
#    $1 -kP 5 -k 5 -D 2e-7 -m 4 -p 2 -s $s -P 1 -C 10000 -o . -c ../../coloring.txt
#done
