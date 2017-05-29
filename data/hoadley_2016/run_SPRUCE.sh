#!/bin/bash
if [ ! $# -eq 2 ]
then
    echo "Usage: $0 <enumerate_executable> <visualize_executable>" >&2
    exit 1
fi

for f in {A1,A7}
do
    echo "Solving patient $f..."
    $1 $f/raw/${f}_spruce_betas.0.9999.tsv > $f/$f.spruce.solution
    python transform.py $f/raw/${f}_spruce_betas.0.9999.tsv > $f/F.tsv
    cd $f
    $2 -t $f.spruce.solution
    cd ..
done
