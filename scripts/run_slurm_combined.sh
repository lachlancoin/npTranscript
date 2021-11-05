#!/bin/bash
#SBATCH --job-name=npTr_comb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16000 # mb
#SBATCH --time=100:00:00
#SBATCH --output=%x.%j.stdout
#SBATCH --error=%x.%j.stderr
#SBATCH --cpus-per-task=8
#SBATCH --partition=physical,snowy


## Load required modules
module load gcc/8.3.0
module load samtools/1.9
module load minimap2/2.21
#module load minimap2/2.17
module load java

##example
# sbatch run_slurm_combined.sh combined todo_fusion.txt opts_run.txt resdir optsType=sc_DNA opts_annot.txt
# sbatch run_slurm_combined.sh combined todo_promethion.txt opts_run_rna.txt --optsType=sc_cDNA
# sbatch --array=2 run_slurm_combined.sh combined todo_promethion.txt opts_run_rna.txt --optsType=sc_cDNA
#find /data/gpfs/projects/punim1068/host_analysis/direct_RNA_sequin_genome/Vero/ -name *bam > todo.txt
export JSA_MEM=12000m
mm2_mem="15g"   
mm2_threads=8
mm2_batch=5M

m=$SLURM_ARRAY_TASK_ID
if [ ! $m ]; then
 m=1
fi


if [ ! $npTranscript ] ; then
export	npTranscript=${HOME}/github/npTranscript
#export npTranscript=/data/gpfs/projects/punim1068/npTranscript_streaming
fi 


species=$1
shift;
if [ ! $species ]; then 
	echo "first option needs to be 'monkey'  or ' human' "
	exit;
fi

db_path=/data/gpfs/projects/punim1068/db
mm2_splicing="uf"  # for host
# mm2_splicing=$(grep 'mm2_splicing' ${opts_host} | grep -v '^#'  | tr -s ' ' | cut -f 1 -d ' ' | cut -f 2 -d '=')
	
if [ $species == "human" ]; then
	reference="${db_path}/merged/sequin_genome/VIC_virus_human_sequin_genome_ensembl_pri_merged.fasta"
	coord_file="${db_path}/merged/human_ensembl_sequin_merged.gtf"
elif [ $species == "human3" ]; then
        reference="./human3.fa.gz"
elif [ $species == "monkey" ]; then
	reference="${db_path}/merged/sequin_genome/VIC_virus_monkey_sequin_genome_merged.fasta"
	coord_file="${db_path}/merged/monkey_sequin_merged.gtf"
elif [ $species == "human_only" ]; then
	reference="./GRCh38_full_analysis_set_plus_decoy_hla.fa.gz"
	coord_file="./Homo_sapiens.GRCh38.100.gtf.gz"
elif [ $species == "SARSCov2" ]; then
	mm2_splicing="un"
	reference="${npTranscript}/data/SARS-Cov2/wuhan/wuhan_coronavirus.fasta.gz"
	coord_file="${npTranscript}/data/SARS-Cov2/wuhan/Coordinates.csv"
elif [ $species == "rRNA" ]; then
	reference="${npTranscript}/data/U13369/Human_ribosomal_DNA_complete_repeating_unit.fasta.gz"
	#reference="${npTranscript}/data/U13369/18S_28S.fa.gz"
elif [ $species == "combined" ]; then
	reference="${npTranscript}/data/U13369/combined.fa.gz"
elif [ $species == "combinedT2T" ]; then
        reference=./combinedT2T.fa.gz
else 
	echo "first option needs to be 'monkey'  or ' human' "
	exit
fi

if [ ! -s ${reference}.mmi ]; then
  echo "${reference}.mmi does not exist" ;
  exit 1;
fi
echo $reference


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
if [ ! $opts_annot ];then
 opts_annot="opts_annot.txt"
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


#resdir="./res_${species}_${todo}.${n}"  #.${n}_${dat}_${SLURM_JOB_ID}"
#mkdir $resdir
#cat $todo | head -n $n | tail -n 1 > ${resdir}/${todo}.$n
#cp ${opts_host}  ${resdir}


mainclass="npTranscript.run.ViralTranscriptAnalysisCmd2"
JSA_CP=${npTranscript}/java/target/npTranscript-1.0.jar:${npTranscript}/java/target/classes/
#str="java -Xmx${JSA_MEM} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -classpath ${JSA_CP} ${mainclass} $@"

echo "resdir=${resdir}"

   fastqf=$(find ${resdir} -name "x*.fastq" | sort | head -n $m | tail -n 1)  
   echo "new fastq ${fastqf}"
    readsOutput=${fastqf}.reads.txt.gz
   str="java -Xmx${JSA_MEM} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -classpath ${JSA_CP} ${mainclass}  --reference=${reference} --readsOutputFile=${readsOutput} --resdir=${resdir} --optsFile=${opts_host} --inputFile=${fastqf}.bam ${@}" 
   echo $str
   echo $@
   if [ ! -s ${fastqf}.bam ] ; then
 #   minimap2 -t ${mm2_threads} -ax splice -${mm2_splicing} -K ${mm2_batch} ${reference}.mmi ${fastqf}  | java -Xmx${JSA_MEM} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -classpath ${JSA_CP} ${mainclass} ${barcode_opt} --reference=${reference} --resdir=${resdir} --readsOutputFile=${readsOutput} --optsFile=${opts_host} - 
     minimap2 -t ${mm2_threads} -ax splice -${mm2_splicing} -K ${mm2_batch} ${reference}.mmi ${fastqf} -o ${fastqf}.bam  
     cat /dev/null > ${fastqf}
     #rm -f ${fastqf}
   fi 
#   $str
   if [ ! -s ${readsOutput} ]; then
    java -Xmx${JSA_MEM} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -classpath ${JSA_CP} ${mainclass}  --reference=${reference} --readsOutputFile=${readsOutput} --resdir=${resdir} --optsFile=${opts_host} --inputFile=${fastqf}.bam $@
    if [ -s ${readsOutput} ]; then 
       rm -f ${fastqf}.bam
    fi
   fi

  #java -Xmx${JSA_MEM} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -classpath ${JSA_CP} ${mainclass} ${barcode_opt} --reference=${reference} --overwrite=false --resdir=${resdir} --optsFile=${opts_host} --inputFile=${fastqf}.bam
#echo " running java" 
# java -Xmx${JSA_MEM} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -classpath ${JSA_CP} ${mainclass} --barcode_file=${barcode} --reference=${reference} --resdir=${resdir} --optsFile=${opts_host} --inputFile=${resdir}/out.bam $@
 # echo $str
 # $str
  #exit 1;

### default not to run this should move it elsewhere

if [ ! -f $opts_annot ]; then
  echo "need to have a file opts_annot.txt for annotation step"
  exit ;
fi
cp ${opts_annot} ${resdir}
cd ${resdir}
h5file="0.isoforms.h5"

if [ ! -f $h5file ]; then
  echo "h5 file not created"
  exit 1;
fi

if [ ! $coord_file ]; then
  echo "no coordinate file"
  exit 1;
fi

if [ $opts_annot ]; then
  bash ${npTranscript}/scripts/run_annot.sh  --reference=${reference} --annotation=${coord_file} --resdir="./" --optsFile=${opts_annot} --h5file=${h5file}
fi



