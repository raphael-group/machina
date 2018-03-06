#!/bin/bash
if [ ! -e sciclone ]
then
    mkdir sciclone
fi

for m in {m5,m8}
do
    if [ ! -e sciclone/$m ]
    then
        mkdir sciclone/$m
    fi
    for p in {mS,S,M,R}
    do
        if [ ! -e sciclone/$m/$p ]
        then
            mkdir sciclone/$m/$p
        fi
        for f in ../../../data/sims/${m}/$p/reads_seed*.tsv
        do
            s=$(basename $f .tsv | sed -e s/reads_seed//g)
            if [ ! -e sciclone/$m/$p/$s ]
            then
                mkdir sciclone/$m/$p/$s
            fi
        done
    done
done

for m in {m5,m8}
do
    for p in {mS,S,M,R}
    do
        for f in ../../../data/sims/$m/$p/reads_seed*.tsv
        do
            s=$(basename $f .tsv | sed -e s/reads_seed//g)
            echo Generating PyClone slurm file for $m, seed $s, pattern $p...
            f="sciclone_slurm/${m}_${p}_${s}.slurm"
            echo "#!/bin/bash" > $f
            echo "#SBATCH --job-name=\"sci_${m}_${p}_${s}\"" >> $f
            echo "#SBATCH -c 1" >> $f
            echo "#SBATCH -t 48:00:00" >> $f
            echo "#SBATCH --mem=2G" >> $f
            echo "cd /n/fs/ragr-data/projects/MACHINA/machina_rev/result/sims/clustering" >> $f
            echo "python ./convert_reads_to_sciclone.py ../../../data/sims/$m/$p/reads_seed${s}.tsv sciclone/$m/$p/$s/" >> $f
            echo "Rscript sciclone.R sciclone/${m}/${p}/$s 10" >> $f
        done
    done
done
