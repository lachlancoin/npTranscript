#!/bin/bash
# n=$(wc -l todo.txt)
# sbatch --array 1 consensus.sh
# sbatch --array 1  ../consensus_ORF.sh  0  ID
# sbatch --array 1  ../consensus_ORF.sh  0  ORF
#clear previous results:
# pwd=$(pwd); while read line; do cd $line; ls | grep results | xargs rm -rf ; cd $pwd; done < todo.txt

#SBATCH --job-name=abpoa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=31000 # mb
#SBATCH --time=100:00:00
export JSA_MEM=8000m
abpoa='/sw/abpoa/v1.0.1/abpoa'
npTranscript=${HOME}/github/npTranscript

chrom=$1
ID=$2
reference=$3
if [ ! $1 ]; then
 chrom="0"
fi


if [ ! $2 ]; then
 ID=".fa"
fi

if [ ! $3 ] ; then
reference="${npTranscript}/data/SARS-Cov2/VIC01/wuhan_coronavirus_australia.fasta.gz"
fi

n=$SLURM_ARRAY_TASK_ID
if [ ! $SLURM_ARRAY_TASK_ID ] ; then
 n=$4
fi

if [ ! $n ]; then
 n=1
fi

echo "reference "$reference

zipfile="${chrom}.clusters.zip"
tododir="${chrom}.clusters"
outdir="${chrom}.consensus"
mkdir -p $tododir
mkdir -p $outdir
todolist=$chrom.$n.$ID.todo.txt
fulllist=$(less $zipfile  | grep $ID  | grep .fa | rev | cut -f 1 -d ' ' | rev)
#echo $fulllist
echo $fulllist | wc
#entry=$(less ${chrom}.clusters.zip  | grep $ID | grep .fa  | rev | cut -f 1 -d ' ' | rev | tail -n +$n | head -n 1)
entry=$(less ${chrom}.todo.txt | rev | cut -f 1 -d ' ' | rev | tail -n +$n | head -n 1)
ID1=$(echo $entry | cut -f 1 | rev | cut -f 2- -d '.' | rev)
outfile="${outdir}/${ID1}.fasta"
	echo $entry
	entryname=$(echo $entry | sed 's/.fa//')
	clustid=$(echo $entry | cut -f 1-2 -d .)
	#desc=$(less $chrom.transcripts.txt.gz | grep -P "${clustid}\t" | head -n 1 | sed 's/\t/;/g' | cut -f 1-12 -d ';')
	tmpfile="${tododir}/${entry}"
	unzip -p $zipfile $entry > $tmpfile
	firstline=$(head -n 1 $tmpfile |  cut -f 2- -d ' ')
	count=$(grep '^[>@]' $tmpfile | wc -l)

#filter to improve consensus  ## sd_thresh tolerance min_seqs
	bash ${npTranscript}/scripts/run_filter.sh   $tmpfile 1.0 5 1

	${abpoa} $tmpfile | sed "s/Consensus_sequence/${entryname} ${firstline} ${count}/" >> $outfile
	#rm -f $tmpfile
#rm $todolist

#rmdir $tododir


##NOW MAP TO REF
export JSA_MEM=8000m
bash ${npTranscript}/scripts/run_map_consensus.sh   --reference=${reference} --fasta=${outfile} --offset 100


