#!/bin/bash
for m in {m5,m8}
do
#m=m8
for p in {mS,S,M,R}
do
    for f in ../../../../data/sims/${m}/$p/reads_seed*.tsv
    do
	s=$(basename $f .tsv | sed -e s/reads_seed//g)
	echo Generating TreeOmics slurm file for seed $s pattern $p...
	
	f="${p}_${s}_${m}.slurm"
	echo "#!/bin/bash" > $f
	echo "#SBATCH --job-name=\"${p}_${s}_treeomics\"" >> $f
	echo "#SBATCH -c 16" >> $f
	echo "#SBATCH -t 1:30:00" >> $f
	echo "#SBATCH --mem=100G" >> $f
        echo "cd /n/fs/ragr-data/projects/MACHINA/treeomics/src" >> $f
        echo "source activate treeomics_env" >> $f
	echo "python treeomics -r ../../machina/result/sims/treeomics/input_${m}/${p}_seed${s}_variant.txt -s ../../machina/result/sims/treeomics/input_${m}/${p}_seed${s}_coverage.txt -o ${m}/output/${p}seed${s} -t 3600 " >> $f

	ff="${p}_${s}_${m}_SUBCLONES.slurm"
	echo "#!/bin/bash" > $ff
	echo "#SBATCH --job-name=\"${p}_${s}_treeomics_SUBCLONES\"" >> $ff
	echo "#SBATCH -c 16" >> $ff
	echo "#SBATCH -t 4:30:00" >> $ff
	echo "#SBATCH --mem=200G" >> $ff
        echo "cd /n/fs/ragr-data/projects/MACHINA/treeomics/src" >> $ff
        echo "source activate treeomics_env" >> $ff
	echo "python treeomics -r ../../machina/result/sims/treeomics/input_${m}/${p}_seed${s}_variant.txt -s ../../machina/result/sims/treeomics/input_${m}/${p}_seed${s}_coverage.txt -o ${m}/outputSUB/${p}_seed${s}SUBCLONES -u -t 10800" >> $ff
    done
done
done
