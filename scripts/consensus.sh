alias abpoa='/sw/abpoa/v1.0.1/abpoa'
chrom=$1
zipfile="${chrom}.clusters.zip"
tododir=$(echo $todo | sed 's/.zip//')
mkdir -p $tododir
todolist=$chrom.todo.txt
unzip -l $zipfile  | tr -s ' ' | cut -f 5 -d ' ' | grep '.fa' > $todolist
outfile="${chrom}.consensus.fasta"
while read entry; do
	echo $entry
	#entry=$(head -n $n $todolist | tail -n +1)
	entryname=$(echo $entry | sed 's/.fa//')
	clustid=$(echo $entry | cut -f 1-2 -d .)
	desc=$(less $chrom.transcripts.txt.gz | grep -P "${clustid}\t" | head -n 1 | sed 's/\t/;/g' | cut -f 1-12 -d ';')
	tmpfile="${tododir}/${entry}"
	#outfile="${tododir}/${entry}.out"
	unzip -p $zipfile $entry > $tmpfile
	abpoa $tmpfile | sed "s/Consensus_sequence/${entryname} ${desc}/" >> $outfile
done < $todolist
rm $todolist
rm -rf $tododir
