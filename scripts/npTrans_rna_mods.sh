#!/bin/bash -e

#conda activate ~/nanoraw_env

###REQUIRES:
#	tombo
#	seqtk
#	SquiggleKit

###INPUTS:
#	raw_reads_basedir $1
#       npTranscript cluster.zip file $2
#	bascalled reads $3
#	statistics-file-basename/experiment name $4
#       consensus_dir $5

consensus_ref=$5
if [ ! $5 ] ; then
	consensus_ref="/DataOnline/Data/raw_external/Coronavirus/rerio_basecalled/double-splice/results_Taiaroa/0.Taiaroa.consensus/wuhan_coronavirus_australia" ;
fi

#store parameters
echo ${@} > $4_params.txt

#list of target ORF clusters
transcripts=($(unzip -l $2 | awk '{print $4}'| grep 'leader;leader,[A-Z]*[0-9]*[a-z]*;end'))
#sort into unique
uniq_ORFs=(`echo ${transcripts[@]} | awk 'BEGIN{FS=";";RS=" "} {print $2}' | sort | uniq`)
echo ${uniq_ORFs[@]}

echo 'seqtk sample'
cat_readfile=$4_subseq.fq
for uo in ${uniq_ORFs[@]}; do target_files=(`grep $uo <(echo ${transcripts[@]} | sed 's/ /\n/g') `) ;
#subseq_readfile=subseq_$uo_"$4".fastq;
#uniq_cluster[$uo]=$subseq_readfile ;
echo ${target_files[@]} >> "$4"_uotargets.txt ;
unzip -p $2 ${target_files[@]} | awk 'BEGIN{RS=">"} {print $1}' | shuf -n 100 - | tee "$4"_sampled_seqs.txt | seqtk subseq $3 - ;
samps=$(wc -l "$4"_sampled_seqs.txt | cut -d " " -f 1); if [ $samps -lt 100 ]
	then printf "$uo\t$samps\n" >> $4_null_samples.txt 
fi ;
done >> $cat_readfile

#fetch matching fast5 files
echo 'fast5 fetcher'
raw_dir="$4"_fast5
mkdir $raw_dir

python3 ~/TEST_squiggle/SquiggleKit/fast5_fetcher_multi.py -q $cat_readfile  \
-s $(dirname $3)/sequencing_summary* -o $raw_dir -m $1

echo 'tombo preprocess'
tombo preprocess annotate_raw_with_fastqs --overwrite --fast5-basedir $raw_dir --fastq-filenames $cat_readfile

echo 'tombo resquiggle'
for transcript in ${uniq_ORFs[@]}; do
cat `find $consensus_ref -name *$transcript";end"*.ref.fa | head -n 1` ; done >> $4_tombo_reference.fa
grep '>' $4_tombo_reference.fa
tombo resquiggle --rna --overwrite --ignore-read-locks $raw_dir $4_tombo_reference.fa

echo 'tombo detect modifications'
tombo detect_modifications alternative_model --rna --fast5-basedirs $raw_dir --statistics-file-basename $4 --alternate-bases 5mC

echo 'tombo text output'
tombo text_output browser_files --browser-file-basename $4 --file-types dampened_fraction  \
--statistics-filename $4.5mC.tombo.stats --fast5-basedirs $raw_dir
