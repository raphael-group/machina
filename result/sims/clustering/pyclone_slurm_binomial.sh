#!/bin/bash
for m in {m5,m8}
do
    for p in {mS,S,M,R}
    do
        for f in ../../../../data/sims/$m/$p/reads_seed*.tsv
        do
            s=$(basename $f .tsv | sed -e s/reads_seed//g)
            echo Generating PyClone slurm file for seed $s pattern $p...
            
            f="${m}_${p}_${s}.slurm"
            echo "#!/bin/bash" > $f
            echo "#SBATCH --job-name=\"bin_${m}_${p}_${s}\"" >> $f
            echo "#SBATCH -c 1" >> $f
            echo "#SBATCH -t 24:00:00" >> $f
            echo "#SBATCH --mem=4G" >> $f
            echo "cd /n/fs/ragr-data/projects/MACHINA/machina_rev/result/sims/clustering/pyclone/$m/$p/$s/" >> $f
            echo "PyClone run_analysis --config_file config_binomial.yaml" >> $f
            echo "PyClone build_table --burnin 1000 --config_file config_binomial.yaml --out_file clustering_binomial.tsv --table_type loci" >> $f
        done
    done
done
