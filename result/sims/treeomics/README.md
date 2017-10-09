To run `Treeomics`:

1. Generate Treeomics input files using `generate_input.sh`

2. Generate Treeomics slurm files:
```
    mkdir slurm
    cd slurm
    ../generate_slurm.sh
```
3. Run Treeomics slurm scripts
```
    for f in *.slurm; do sbatch $f; done
```
4. Process Treeomics results
```
    ./process_m5_noSub.sh ../../../build/visualizeclonetree ../../../build/rf
    ./process_m5_sub.sh ../../../build/visualizeclonetree ../../../build/rf
    ./process_m8_noSub.sh ../../../build/visualizeclonetree ../../../build/rf
    ./process_m8_sub.sh ../../../build/visualizeclonetree ../../../build/rf
```
