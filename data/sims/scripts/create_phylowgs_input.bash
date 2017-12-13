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
        dir="../${m}/${p}/phylowgs"
        if [ ! -d $dir ]; then
            echo "Creating directory: $dir"
            mkdir $dir
        else
            echo "Directory exists: $dir"
        fi
        
        for filename in `ls ../${m}/${p}/reads*`;
        do
            echo "Reading: ${filename}"
            python create_phylowgs_input.py $filename $dir
        done
        
    done
done 
