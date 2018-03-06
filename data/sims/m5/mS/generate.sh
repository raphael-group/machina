#!/bin/bash
# run on turing.cs.princeton.edu

if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <simulate_executable>" >&2
    exit 1
fi

#$1 -kP 2 -k 1 -D 2e-7 -m 4 -p 0 -s 0 -N 10 -C 200 -o . -c ../../coloring.txt

for C in {60,200,10000}
do
    if [ ! -e C${C} ]
    then
        mkdir C${C}
    fi
    for P in {0.6,0.8,1}
    do
        if [ ! -e C${C}/P${P} ]
        then
            mkdir C${C}/P${P}
        fi
	for k in {1,3,5}
	do
            if [ ! -e C${C}/P${P}/k${k} ]
            then
                mkdir C${C}/P${P}/k${k}
            fi
            for s in {0,10,12,2,3,4,5,7,8,9}
            do
                if [ $k -eq 1 ]
                then
                    $1 -kP 2 -k 1 -D 2e-7 -m 4 -p 0 -s $s -C $C -E 0.001 -P $P -o C${C}/P${P}/k${k} -c ../../coloring.txt
		else
                    $1 -kP $k -k $k -D 2e-7 -m 4 -p 0 -s $s -C $C -E 0.001 -P $P -o C${C}/P${P}/k${k} -c ../../coloring.txt
                fi
            done
	done
    done
done
