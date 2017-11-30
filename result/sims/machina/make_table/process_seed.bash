## Get seeds
seed=$2
pattern=$1
nsites=$3

if [ -z "$pattern" ]
then
    echo "\\hline "
    exit 0
fi
#echo $seed, $pattern, $nsites

dataloc=../../../../data/sims/m${nsites}/${pattern}/
resultloc=../
true_graph=${dataloc}/G_seed${seed}.tree

echo $seed
echo " & "

# Get number of sites
cat $true_graph | tr ' ' '\n' | sort | uniq | wc -l 
echo " & "

# get number of mut trees
mut_trees=${resultloc}/mut_trees_m${nsites}/mut_trees_${pattern}_seed${seed}.txt
head -2 $mut_trees | tail -1 | sed 's/#trees//'
echo " & "

if [ $pattern == "mS" ];
then
    echo $pattern
else
    echo p$pattern
fi
echo " & "

# Get number of comgs
echo "("
cat $true_graph | wc -l 
echo ", "
cat $true_graph | sort | uniq | wc -l 
echo ") & "

# Get minimum achieved values 
result_file=${resultloc}/output_m${nsites}/${pattern}_${seed}.txt
python get_min_values.py $result_file
echo " \\\ "
