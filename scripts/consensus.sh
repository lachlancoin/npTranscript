#!/bin/bash
# n=$(wc -l todo.txt)
# sbatch --array 0 consensus.sh

#clear previous results:
# pwd=$(pwd); while read line; do cd $line; ls | grep results | xargs rm -rf ; cd $pwd; done < todo.txt

#SBATCH --job-name=abpoa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=9000 # mb
#SBATCH --time=100:00:00
#SBATCH --output=abpoa_%A_%a.out
#SBATCH --error=abpoa_%A_%a.err

abpoa='/sw/abpoa/v1.0.1/abpoa'
npTranscript=${HOME}/github/npTranscript
reference="${npTranscript}/data/SARS-Cov2/wuhan_coronavirus_australia.fasta.gz"

chrom=0
n=$SLURM_ARRAY_TASK_ID
#chrom=$SLURM_ARRAY_TASK_ID
#if [ ! $SLURM_ARRAY_TASK_ID ] ; then
#chrom=$1
#fi


zipfile="${chrom}.clusters.zip"
tododir="${chrom}.clusters"
mkdir -p $tododir
todolist=$chrom.$n.todo.txt
unzip -l $zipfile  | tr -s ' ' | cut -f 5 -d ' ' | grep -v 'ID' | grep '.fa' | tail -n +$n | head -n 1 > $todolist
outfile="${chrom}.$n.consensus_orf.fasta"
rm $outfile
while read entry; do
	echo $entry
	#entry=$(head -n $n $todolist | tail -n +1)
	entryname=$(echo $entry | sed 's/.fa//')
	clustid=$(echo $entry | cut -f 1-2 -d .)
	#desc=$(less $chrom.transcripts.txt.gz | grep -P "${clustid}\t" | head -n 1 | sed 's/\t/;/g' | cut -f 1-12 -d ';')
	tmpfile="${tododir}/${entry}"
	#outfile="${tododir}/${entry}.out"
	unzip -p $zipfile $entry > $tmpfile
	${abpoa} $tmpfile | sed "s/Consensus_sequence/${entryname}/" >> $outfile
	rm -f $tmpfile
done < $todolist
rm $todolist

#rmdir $tododir


##NOW MAP TO REF
export JSA_MEM=8000m
bash ${npTranscript}/scripts/run_map_consensus.sh   --reference=${reference} --fasta=${outfile} --offset 10000


