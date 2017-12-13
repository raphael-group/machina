#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <m>" >&2
    exit 1
fi

m=$1

for p in {mS,S,M,R}
do
    for f in ../../../../data/sims/${m}/$p/reads_seed*.tsv
    do
        s=$(basename $f .tsv | sed -e s/reads_seed//g)
        echo Generating MACHINA PMH-TI slurm file for seed $s pattern $p...
        
        f="MACHINA_${m}_${p}_${s}.slurm"
        echo "#!/bin/bash" > $f
        echo "#SBATCH --job-name=\"${p}_${s}_MACHINA_${m}\"" >> $f
        echo "#SBATCH -c 2" >> $f
        echo "#SBATCH -t 1:30:00" >> $f
        echo "#SBATCH --mem=10G" >> $f
        echo "module load boost" >> $f
        echo "module load gurobi" >> $f
        echo "cd /n/fs/ragr-data/projects/MACHINA/machina_rev/result/sims/machina/" >> $f
        echo "mkdir output_${m}/${p}_${s}" >> $f
        echo "../../../build/generatemutationtrees -tl 10 input_${m}/cluster_${p}_seed${s}.tsv > mut_trees_${m}/mut_trees_${p}_seed${s}.txt" >> $f
        echo "../../../build/pmh_ti -p P -c ../../../data/sims/coloring.txt -m 3 -t 2 -o output_${m}/${p}_${s}/ -F input_${m}/cluster_${p}_seed${s}.tsv -barT mut_trees_${m}/mut_trees_${p}_seed${s}.txt > output_${m}/${p}_${s}.txt" >> $f
    done
done
