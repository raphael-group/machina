set -euo pipefail
v=$1

basev=`echo $v | sed "s/..\/..\/..\/..\/data\/sims\///"`
dir=../`dirname $basev`
dir=`echo $dir | sed "s/phylowgs//"`
file=`basename $basev`

pushd $dir
if [ ! -d $file ]; then
    mkdir $file
fi

cd $file

ls *
rm *

touch EMPTY
echo "cnv    a   d   ssms    physical_cnvs" > EMPTY

echo "Running on $v"

phylowgs='../../../../../../../other_software/phylowgs/evolve.py'
python2 $phylowgs ../../$v EMPTY --random-seed 1

#python2 $phylowgs ../../../scripts/my_data.txt EMPTY
#python2 $phylowgs ../../../scripts/their_data.txt EMPTY

popd

