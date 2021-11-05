#!/bin/bash
#SBATCH --job-name=submit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=7000 # mb
#SBATCH --time=10:00
#SBATCH --output=%x.%j.stdout
#SBATCH --error=%x.%j.stderr
#SBATCH --cpus-per-task=8
#SBATCH --partition=physical,snowy


#sbatch submit.sh 1 todo_reference.txt RNA 
#bash submit.sh 1 todo_reference.txt RNA
#bash submit.sh 1 todo_reference.txt RNA
#sbatch submit.sh 4 todo_promethion.txt sc_cDNA
n=$1
if [ ! $n ] ; then 
echo " problem " ;
exit 1;
fi

shift;


batch_size=1000000


species="combinedT2T"
todo="todo_promethion.txt"
opts_host="opts_run_rna.txt"
optsType="sc_cDNA"

if [ $1 ]; then
 todo=$1;
 shift;
fi

if [ $1 ]; then
 optsType=$1;
fi

resdir="./res_${species}_${todo}.${n}" 
echo $resdir
mkdir -p $resdir
fastq=$(cat $todo | head -n $n | tail -n 1)
fqname=$(echo $fastq | rev | cut -f 1 -d '/' | rev)
echo $fqname

cp ${opts_host}  ${resdir}
cat $todo | head -n $n | tail -n 1 > ${resdir}/${todo}.$n

mcomplete=$(find ${resdir} -name "x*.reads.txt.gz" -type f | wc -l )
mbam=$(find ${resdir} -name "x*.bam" -type f | wc -l )

if [ $mcomplete -gt 0 ]; then
 echo " already done  ${mcomplete}" 
	exit 1;
fi


fastq=$(cat $todo | head -n $n | tail -n 1)
  fqname=$(echo $fastq | rev | cut -f 1 -d '/' | rev)
  echo $fqname
  barcode_list=$(grep 'barcode_list' ${opts_host} | grep -v '^#'  | tr -s ' ' | cut -f 1 -d ' ' | cut -f 2 -d '=')
  echo $barcode_list
  #barcode=${barcode_list}
  if [ -s $barcode_list ]; then
#        barcode=$(cat ${barcode_list.txt}  | tr -s ' ' | sed 's/ /     /g' | grep ${fqname} | cut -f 2)
        barcode=$(cat $barcode_list | grep ${fqname} | tr -s ' '  | sed 's/ /	/g'  | cut -f 2)
        if [ $(echo $barcode | wc -c) -gt 1 ]; then
                barcode_opt="--barcode_file=${barcode}"
        else
                barcode_opt=""
        fi
        echo "HH  barcode ${barcode_opt}"
  else
        barcode="none"
        barcode_opt="--barcode_file=none" 
        echo "exiting, no barcode"
        exit 1;
  fi

echo "here 1 barcode_opt  ${barcode_opt}" 

#if [ $2 ]; then
#if [ $2 == "clean" ]; then 
#  echo " cleaning" 
#  find ${resdir} -name "x*.fastq" -type f | xargs rm -f
#  find ${resdir} -name "x*.fastq.bam" -type f | xargs rm -f
#  exit 1;
#fi
#fi
is_ref=0
m=$(find ${resdir} -name "x*.fastq" -type f | wc -l )
if [ $m -eq 0 ]; then
        cd ${resdir}
        if [ $(echo $fastq | grep '.gz' | wc -l ) -eq 0 ]; then
	 if [ $is_ref -eq 0 ] ;then
          cat ${fastq}   | split -l ${batch_size} -d  --additional-suffix=.fastq
	 else
          cat ${fastq} | sed 's/ /_/g'  | split -l ${batch_size} -d  --additional-suffix=.fastq
	 fi
        else
	 if [ $is_ref -eq 0 ] ;then
          	zcat ${fastq}   | split -l ${batch_size} -d  --additional-suffix=.fastq
	 else 
         	zcat ${fastq} | sed 's/ /_/g'   | split -l ${batch_size} -d  --additional-suffix=.fastq
         fi
        fi
#--verbose | cut -f 3 -d ' ' | sed "s/’//g" | sed "s/‘//g" | xargs gzip
        cd ..
	m=$(find ${resdir} -name "x*.fastq" -type f | wc -l )
fi


if [ $3 ]; then
 m=$3
fi

str="sbatch --array=1-${m} run_slurm_combined.sh ${species} ${todo} ${opts_host} ${resdir} --optsType=${optsType}  ${barcode_opt}" 

echo $str
#$str
#if [ $2 ] ;then 
#	$str
#else
	echo $str  >> to_submit.${n}.sh
#fi
