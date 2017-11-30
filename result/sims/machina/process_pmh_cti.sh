#!/bin/bash
if [ ! $# -eq 3 ]
then
    echo "Usage: $0 <RESULT_DIR> <rf_executable> <m>" >&2
    exit 1
fi
m=$3

echo "pattern,seed,mut_tree,enforced,inferred,mu,gamma,sigma,method,RF,recallT,precisionT,FscoreT,recallG,precisionG,FscoreG,recallMultiG,precisionMultiG,FscoreMultiG"
for f in $1/*.txt
do
    d=`basename $f .txt`
    p=$(echo $d | sed -e s/_.*//g)
    s=$(echo $d | sed -e s/.*_//g)
    for t in `python extract_minimum.py $f`
    do
        pattern=$(echo $t | cut -d',' -f2)
        idx=$(echo $t | cut -d',' -f1)
	rf=$(echo -n $($2 ../../../data/sims/${m}/$p/T_seed$s.tree ../../../data/sims/$m/$p/T_seed$s.labeling $1/$d/${idx}-T-P-${pattern}.tree $1/$d/${idx}-T-P-${pattern}.labeling | tail -n 1 | cut -d' ' -f3))

        if [ $p == "mS" ];
        then
            echo -n $p,
        else
            echo -n p$p,
        fi

        echo -n $s,$t,MACHINA,$rf,
        python evaluate_mig_history.py ../../../data/sims/${m}/$p/T_seed${s}.tree ../../../data/sims/$m/$p/T_seed${s}.vertex.labeling ../../../data/sims/$m/$p/G_seed${s}.tree $1/${d}/${idx}-T-P-${pattern}.tree $1/${d}/${idx}-T-P-${pattern}.labeling $1/${d}/${idx}-G-P-${pattern}.tree
       	#tail -n 1 ${p}_seed${s}.RF.txt | cut -d' ' -f3  >> results_m7.txt
    done
done

