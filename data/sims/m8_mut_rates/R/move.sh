#!/bin/bash

for s in {10157,1070,30342,30718,30527,30604,394,9,905,981}
do
    for mut in {0.1,1,10}
    do
        mv ${mut}/G_seed${s}.dot G_seed${s}_mut${mut}.dot
        mv ${mut}/G_seed${s}.tree G_seed${s}_mut${mut}.tree
        mv ${mut}/T_seed${s}.dot T_seed${s}_mut${mut}.dot
        mv ${mut}/T_seed${s}.tree T_seed${s}_mut${mut}.tree
        mv ${mut}/T_seed${s}.labeling T_seed${s}_mut${mut}.labeling
        mv ${mut}/T_seed${s}.vertex.labeling T_seed${s}_mut${mut}.vertex.labeling
        mv ${mut}/T_all_seed${s}.dot T_all_seed${s}_mut${mut}.dot
        mv ${mut}/clustering_observed_seed${s}.txt clustering_observed_seed${s}_mut${mut}.txt
        mv ${mut}/drivers_seed${s}.txt drivers_seed${s}_mut${mut}.txt
        mv ${mut}/reads_seed${s}.tsv reads_seed${s}_mut${mut}.tsv
    done 
done
