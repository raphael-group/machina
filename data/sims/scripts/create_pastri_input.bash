####
#   Converts all input files to phylowgs format
####

Ms='m5 m8'
Ps='M mS R S'

for m in $Ms ;
do
    for p in $Ps ;
    do
        echo "-----------"
        echo "Processing ${m} with pattern ${p}"
        dir="../${m}/${p}/pastri"
        if [ ! -d $dir ]; then
            echo "Creating directory: $dir"
            mkdir $dir
        else
            echo "Directory exists: $dir"
        fi
        
        for filename in `ls ../${m}/${p}/reads*`;
        do
            cluster_file="input_${m}/cluster_${p}_`expr match "$filename" '.*\(seed[0-9]\+\).tsv'`.txt"
            ls $cluster_file
            python create_pastri_input.py $filename $dir $cluster_file
        done
    done
done 
