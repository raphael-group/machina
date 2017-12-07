

filelist='sim_list.txt'

run_phylowgs(){
    #bash run_phylowgs.bash $1 &> `python id.py $1`
    bash run_phylowgs.bash $1 &> logs/`python id.py $1`

}
export -f run_phylowgs

parallel -j 20  --joblog joblog.txt run_phylowgs ::: `cat $filelist` &
