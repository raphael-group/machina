#set e
# pipefail


####
#   Converts all input files to phylowgs format
####

Ms='m5 m8'
Ps='M mS R S'

for m in $Ms ;
do
    for p in $Ps ;
    do

        for filedir in ../${m}/${p}/reads_seed*.phylowgs; do

            echo "-------------------------------"
            echo "Processing: $filedir"

            id=`expr match "$filedir" '.*seed\([0-9]\+\).*'`

            echo $?
            echo $id
            input_file=../../../../data/sims/${m}/${p}/reads_seed${id}.tsv 

			###
			# For each result, we have to run the script to process the result 
			# unzip the file
			# run parse_trees
			
		    	
			write_results="python ../../../../../other_software/phylowgs/write_results.py"
			
			if [ ! -d $filedir/results ]; then
			    mkdir $filedir/results
			
			    mut_file=$filedir/results/result.muts.json
			    result_file=$filedir/results/result.summ.json
			    tree_file=$filedir/results/result.mutass
			    echo "Processing PhyloWGS results"
			    $write_results result $filedir/trees.zip ${result_file}.gz ${mut_file}.gz ${tree_file}.zip
			    echo "Extracting result files"
			
			    gunzip ${result_file}.gz
			    unzip  ${tree_file}.zip -d $filedir/results > /dev/null
			
			    edgelist_file=$filedir/edgelist.txt
			    leaflabel_file=$filedir/leaflabels.txt
			
			    echo "Parsing trees"
			
			    python parse_trees.py $result_file $input_file $edgelist_file $leaflabel_file
            else
                echo "Results exist"

			fi
        done
    done
done 

