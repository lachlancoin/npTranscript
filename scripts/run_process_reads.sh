#!/bin/bash
#SBATCH --job-name=npTr_proc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32000 # mb
#SBATCH --time=3:59:00
#SBATCH --output=%x.%j.stdout
#SBATCH --error=%x.%j.stderr
#SBATCH --cpus-per-task=8
#SBATCH --partition=physical,snowy


## Load required modules
module load java

##example
#sbatch run_process_reads.sh res_combinedT2T_todo_reference.txt.1 res_combinedT2T_todo_promethion.txt.2 --maxcells=3.049e6 --num_barcodes=3049 --type=spliced
##type is some combination of chrom;nostrand;truncated;nostrand_trunc;spliced;all 
#outf="out.h5"
export JSA_MEM=26000m

suffix=".reads.txt.gz"

if [ ! $npTranscript ] ; then
export	npTranscript=${HOME}/github/npTranscript
#export npTranscript=/data/gpfs/projects/punim1068/npTranscript_streaming
fi 


dat=$(date +%Y%m%d%H%M%S)



if [ ! $(echo $1 | grep '\-\-') ]; then
 refdir=$1
 resdir=$2
 shift;
 shift;
fi



if [ ! $resdir ]; then
 echo " need resdir" 
 exit 1;
fi

if [ ! $refdir ]; then
 echo " need refdir"
 exit 1;
fi


echo "resdir=${resdir}"
echo "refdir=${refdir}"
#cd $resdir

mainclass="npTranscript.run.ProcessReadFile"
JSA_CP=${npTranscript}/java/target/npTranscript-1.0.jar:${npTranscript}/java/target/classes/
  
 java -Xmx${JSA_MEM} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -classpath ${JSA_CP} ${mainclass} --refFileDir=${refdir} --inputFileDir=${resdir} --suffix=${suffix}  $@
