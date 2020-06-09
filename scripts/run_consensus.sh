# this is run in the results director from npTranscript
npTranscript=${HOME}/github/npTranscript
reference=$1
thresh=20
if [ ! $1 ]; then
 reference="${npTranscript}/data/SARS-Cov2/VIC01/wuhan_coronavirus_australia.fasta.gz:${npTranscript}/data/SARS-Cov2/wuhan/wuhan_coronavirus.fasta.gz"
fi

fa=$2
if [ ! $2 ]; then
 fa="clust"
fi


todo=$(ls *.clusters.zip | grep $fa | xargs -I {} echo {} | rev | cut -f 3- -d . | rev)
for i in $todo; do
	echo $i.clusters.zip
	echo $i.clusters
	mkdir $i.clusters
	unzip -n $i.clusters.zip -d $i.clusters
	cd $i.clusters
	for j in *.fa; do
		a=$(grep '>' $j | wc -l)
		if [ $a -lt $thresh ] ;then rm $j ; fi;
	done
	ls | grep '.fa' > ../$i.todo.txt
        cd ../
	#less $i.clusters.zip | grep $fa > $i.todo.txt
	cnt=$(wc -l $i.todo.txt | cut -f 1 -d ' ')
	#cnt=$(less $i.clusters.zip | grep $fa  | wc -l)
  	echo $i $cnt
  	#ls $i.clusters.zip
 	sbatch --array 1-${cnt} ${npTranscript}/scripts/consensus.sh $i $reference
done
