#!/bin/bash
# run on turing.cs.princeton.edu

if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <simulate_executable>" >&2
    exit 1
fi

$1 -k 1 -kP 2 -C 200 -D 2e-7 -m 7 -p 3 -s 0 -N 10 -o . -c ../../coloring.txt

#for s in {109,117,165,22,43,5,61,69,74,9}
#do
#	$1 -p 0 -s $s -o . -c ../coloring.txt
#done
