#!/bin/bash
if [ ! -e clomial ]
then
    mkdir clomial
fi

for m in {m5,m8}
do
    if [ ! -e clomial/$m ]
    then
        mkdir clomial/$m
    fi
    for p in {mS,S,M,R}
    do
        if [ ! -e clomial/$m/$p ]
        then
            mkdir clomial/$m/$p
        fi
        for f in ../../../data/sims/${m}/$p/reads_seed*.tsv
        do
            s=$(basename $f .tsv | sed -e s/reads_seed//g)
            if [ ! -e clomial/$m/$p/$s ]
            then
                mkdir clomial/$m/$p/$s
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
            echo Generating Clomial slurm file for $m, seed $s, pattern $p...
            f="clomial_slurm/${m}_${p}_${s}.slurm"
            echo "#!/bin/bash" > $f
            echo "#SBATCH --job-name=\"clm_${m}_${p}_${s}\"" >> $f
            echo "#SBATCH -c 1" >> $f
            echo "#SBATCH -t 15:00:00" >> $f
            echo "#SBATCH --mem=2G" >> $f
            echo "cd /n/fs/ragr-data/projects/MACHINA/machina_rev/result/sims/clustering" >> $f
            echo "python ./convert_reads_to_clomial.py ../../../data/sims/$m/$p/reads_seed${s}.tsv clomial/$m/$p/$s/" >> $f
            echo "echo \"clomial/$m/$p/$s/\"" >> $f
            echo "Rscript clomial.R clomial/${m}/${p}/$s 10" >> $f
        done
    done
done
