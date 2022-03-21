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


##CONVERT TO MATRIZ
#module load gcc/8.3.0
#module load openmpi/3.1.4
#module load r-bundle-bioconductor/3.9-r-3.6.0
#cd $resdir
#Rscript ${npTranscript}/R/convertHDF5Matrix.R


## Load required modules
module load java
module load gcc/8.3.0
##example
#sbatch run_process_reads.sh res_combinedT2T_todo_promethion.txt.2 --maxcells=3.049e6 --num_barcodes=3049 --type=spliced
#run within resdir
#sbatch run_process_reads.sh resdir_new --maxcells=3.049e6 --num_barcodes=5000     --type=chrom
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
# refdir=$1
 resdir=$1
 shift;
# shift;
fi



if [ ! $resdir ]; then
 echo " need resdir" 
 exit 1;
fi

#if [ ! $refdir ]; then
# echo " need refdir"
# exit 1;
#fi

#mkdir $resdir

#find . -name '*reads.txt.gz' > todo.txt
#while read line; do
#  line1=$(echo $line | sed -s 's/.gz//g')
#  outf1=${resdir}/${line1}
#echo $outf1
#  if [ ! -s $outf1 ]; then
#    echo "zcat ${line} | head -n -5 > ${outf1}"
#    zcat $line | wc -l 
#    zcat $line | head -n -5 > $outf1
#    echo "gzip ${outf1}"
#    gzip $outf1
#  fi
#done < todo.txt



echo "resdir=${resdir}"
#echo "refdir=${refdir}"
#cd $resdir

mainclass="npTranscript.run.ProcessReadFile"
JSA_CP=${npTranscript}/java/target/npTranscript-1.0.jar:${npTranscript}/java/target/classes/
#--inputFileDir=resdir1   --overwrite=true --maxcells=1e5 --num_barcodes=5000 --outputFile=out1.h5    --suffix=.reads.txt.gz --type=chrom

  
# java -Xmx${JSA_MEM} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -classpath ${JSA_CP} ${mainclass} --refFileDir=${refdir} --inputFileDir=${resdir} --suffix=${suffix}  $@
 java -Xmx${JSA_MEM} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -classpath ${JSA_CP} ${mainclass} --inputFileDir=${resdir} --suffix=${suffix}  $@


