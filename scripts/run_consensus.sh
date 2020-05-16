# this is run in the results director from npTranscript

npTranscript=${HOME}/github/npTranscript
cnt=$(less 0.clusters.zip | grep '.fa'  | grep 'ID' | wc -l)
sbatch --array 1-${cnt} ${npTranscript}/scripts/consensus.sh 0 ID

cnt=$(less 0.clusters.zip | grep '.fa'  | grep 'ORF' | wc -l)
sbatch --array 1-${cnt} ${npTranscript}/scripts/consensus.sh 0 ORF

