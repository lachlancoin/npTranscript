##script for running 
##tip - use symbolic link to put this in the directory with bam files
npTranscript=${HOME}/github/npTranscript

bamdir="."
bamfiles=$(ls ${bamdir} | grep '.bam' | xargs -I {} echo ${bamdir}/{})
bamfiles1=$(echo $bamfiles | sed 's/ /:/g')

reference="${npTranscript}/data/SARS-Cov2/wuhan_coronavirus_australia.fasta.gz"
coord_file="${npTranscript}/data/SARS-Cov2/Coordinates.csv"

dat=$(date +%Y%m%d%H%M%S)
resdir="results_${dat}"

#export RESULTS_DIR=$resdir

opts="--bin 10 --breakThresh 1000"
bash ${npTranscript}/scripts/run.sh --bamFile=${bamfiles1}   --reference=${reference}   --annotation ${coord_file}   --resdir ${resdir} ${opts}

