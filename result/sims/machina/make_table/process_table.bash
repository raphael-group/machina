
file=$1

cat header.txt

while read line
do
    #echo $line 
    bash process_seed.bash $line | tr -d '\n' 
    echo
done < $file

cat footer.txt
