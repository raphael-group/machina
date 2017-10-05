#!/bin/bash
for f in raw/A*.txt
do
    echo Processing $f...
    python generateFrequencies.py raw/samples.txt $f prostate 0.025 > `basename $f .txt`.tsv
done
