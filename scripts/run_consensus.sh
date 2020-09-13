# this is run in the results director from npTranscript
npTranscript=${HOME}/github/npTranscript
reference=$2
thresh=20
if [ ! $2 ]; then
 reference="${npTranscript}/data/SARS-Cov2/VIC01/wuhan_coronavirus_australia.fasta.gz:${npTranscript}/data/SARS-Cov2/wuhan/wuhan_coronavirus.fasta.gz"
fi

fa=$1
if [ ! $1 ]; then
 fa="clust"
fi

max_clust=$3
if [ ! $3 ]; then
  max_clust=200
fi

unzip=true

todo=$(ls *.clusters.zip | grep $fa | xargs -I {} echo {} | rev | cut -f 3- -d . | rev)
for i in $todo; do
	echo $i.clusters.zip
	echo $i.clusters
	if [ "$unzip" == true ]; then
		mkdir $i.clusters
		unzip -n $i.clusters.zip -d $i.clusters
	fi;

	cd $i.clusters
	for j in *.fa; do
		a=$(grep '>' $j | wc -l)
		if [ $a -lt $thresh ] ;then rm $j  
		else echo ${j} ${a}
	fi;
	done | sort -rn -k2 | head -n ${max_clust} | cut -d " " -f 1 > ../$i.todo.txt


	#ls | grep '.fa' > ../$i.todo.txt
        cd ../
	#less $i.clusters.zip | grep $fa > $i.todo.txt
	cnt=$(wc -l $i.todo.txt | cut -f 1 -d ' ')
	#cnt=$(less $i.clusters.zip | grep $fa  | wc -l)
  	echo $i $cnt
  	#ls $i.clusters.zip
 	sbatch --array 1-${cnt} ${npTranscript}/scripts/consensus.sh $i $reference
done
