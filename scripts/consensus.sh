#!/bin/bash
#SBATCH --job-name=abpoa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=31000 # mb
#SBATCH --time=100:00:00

#SBATCH --error=%x.%j.stderr
#SBATCH --cpus-per-task=8
#SBATCH --partition=physical,snowy


## Load required modules
#module load anaconda3/5.3.0
#source activate abpoa
#abpoa=abpoa
abpoa=/home/lachlan/github/abPOA-v1.2.5/bin/abpoa

dir=$1.tmp
unzip $1 -d $dir
outf=$1.consensus.fa
rm $outf
todo=$1.todo.txt
find $dir -name '*.fa' > $todo
while read file; do
 echo $file
 filen=$(basename $file)
 echo $filen
 $abpoa $file  | sed "s/Consensus_sequence/${filen}/g" >> $outf 
done < $todo

rm -f $dir/*fa
rmdir $dir
rm $todo
