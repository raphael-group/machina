#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <downsample_executable>" >&2
    exit 1
fi

for s in {1275,1381,1394,1526,1593,1708,1930,219,414,759}
do
    for C in {200,500,1000,10000}
    do 
        for P in {0.5,0.8,1.0}
        do
            for k in {1,2,3}
            do
                $1 reads_seed${s}.tsv -C $C -k $k -P $P -E 0.001 > reads_seed${s}_k${k}_P${P}_C${C}.tsv
            done
        done
    done
done

