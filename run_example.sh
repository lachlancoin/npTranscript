##script for running 
##tip - use symbolic link to put this in the directory with bam files
npTranscript=${HOME}/github/npTranscript

bamfiles=bins/allreads_fq.bam:virion/sorted.virion_refmap.bam
reference="${npTranscript}/data/SARS-Cov2/wuhan_coronavirus_australia.fasta.gz"
coord_file="${npTranscript}/data/SARS-Cov2/Coordinates.csv"


resdir="results_1"
opts="--bin 10 --breakThresh 1000 --cluster_by_annotation true"
bash ${npTranscript}/run.sh --bamFile=${bamfiles}   --reference=${reference} --maxReads 100   --annotation ${coord_file}   --resdir ${resdir} ${opts}

bash ${npTranscript}/run_extract_cluster.sh --inDir ${resdir} 
