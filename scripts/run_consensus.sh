# this is run in the results director from npTranscript
todo=$(ls *.clusters.zip | xargs -I {} echo {} | rev | cut -f 3- -d . | rev)
npTranscript=${HOME}/github/npTranscript
for i in $todo; do
 cnt=$(less $i.clusters.zip | grep '.fa'  | grep 'ORF' | wc -l)
 sbatch --array 1-${cnt} ${npTranscript}/scripts/consensus.sh $i ORF
 cnt=$(less $i.clusters.zip | grep '.fa'  | grep 'ID' | wc -l)
 sbatch --array 1-${cnt} ${npTranscript}/scripts/consensus.sh $i ID
done
