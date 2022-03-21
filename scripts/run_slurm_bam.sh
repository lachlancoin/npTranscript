#!/bin/bash
#SBATCH --job-name=npTr_comb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32000 # mb
#SBATCH --time=100:00:00
#SBATCH --output=%x.%j.stdout
#SBATCH --error=%x.%j.stderr
#SBATCH --cpus-per-task=8
#SBATCH --partition=physical,snowy


## Load required modules
module load gcc/8.3.0
#module load samtools/1.9
#module load minimap2/2.21
#module load minimap2/2.17
module load java

##example
# sbatch run_slurm_combined.sh  todo_promethion.txt opts_run_rna.txt --optsType=sc_cDNA
# sbatch --array=2 run_slurm_combined.sh combined todo_promethion.txt opts_run_rna.txt --optsType=sc_cDNA
export JSA_MEM=12g


m=$SLURM_ARRAY_TASK_ID
if [ ! $m ]; then
 m=1
fi


if [ ! $npTranscript ] ; then
export	npTranscript=${HOME}/github/npTranscript
#export npTranscript=/data/gpfs/projects/punim1068/npTranscript_streaming
fi 


dat=$(date +%Y%m%d%H%M%S)



todo=$1
shift

if [ ! $(echo $1 | grep '\-\-') ]; then
 opts_host=$1
 shift
fi

if [ ! $(echo $1 | grep '\-\-') ]; then
 resdir=$1
 shift
fi


if [ ! $resdir ]; then
 echo " need resdir" 
 exit 1;
fi


if [ ! $todo ];then
 todo="todo.txt"
fi
if [ ! $opts_host ];then
 opts_host="opts_run.txt"
fi

if [ ! -f $opts_host ]; then
  echo "need to have a file opts_run.txt"
  exit ;
fi

if [ ! -f $todo ]; then
  echo "need to have a file todo.txt which lists input files"
  exit ;
fi

dat=$(date +%Y%m%d%H%M%S)



mainclass="npTranscript.run.ViralTranscriptAnalysisCmd2"
JSA_CP=${npTranscript}/java/target/npTranscript-1.0.jar:${npTranscript}/java/target/classes/

echo "resdir=${resdir}"
##
## --readsOutputFile=${readsOutput}
   bamf=$(cat $todo | head -n $m | tail -n 1)
readsO=$(echo $bamf | rev | cut -f 1 -d '/' | rev)
readsOutput="${resdir}/${readsO}.reads.txt.gz"
barcodedir=$(grep ${bamf} barcode_list2.txt | cut -f 2 -d ' ')
echo $barcodedir
java -Xmx${JSA_MEM} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -classpath ${JSA_CP} ${mainclass} --resdir=${resdir} --optsFile=${opts_host} --inputFile=${bamf} --barcode_file=${barcodedir}  --readsOutputFile=${readsOutput} $@





