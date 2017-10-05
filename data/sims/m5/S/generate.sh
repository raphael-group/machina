#!/bin/bash
# run on turing.cs.princeton.edu

if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <simulate_executable>" >&2
    exit 1
fi

$1 -D 2e-7 -k 1 -kP 2 -m 4 -p 1 -s 0 -N 10 -C 200 -o . -c ../../coloring.txt
