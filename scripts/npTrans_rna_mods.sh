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
	consensus_ref="/DataOnline/Data/raw_external/Coronavirus/rerio_basecalled/double-splice/results_Taiaroa/0.0.consensus/wuhan_coronavirus_australia" ;
fi

#store parameters
echo ${@} > $4_params.txt

#make directories. raw_dir for final, resquiggled fast5s, tmp_dir for fast5s being worked on
raw_dir="$4"_fast5
mkdir $raw_dir
tmp_dir=$(mktemp -d)
echo $tmp_dir

#list of target ORF clusters
transcripts=($(unzip -l $2 | awk '{print $4}'| grep 'leader;leader,[A-Z]*[0-9]*[a-z]*;end'))
#sort into unique
uniq_ORFs=(`echo ${transcripts[@]} | awk 'BEGIN{FS=";";RS=" "} {print $2}' | sort | uniq`)
echo ${uniq_ORFs[@]}

cat_readfile=$4_subseq.fq

#START LOOP
for uo in ${uniq_ORFs[@]}; do echo $uo;
echo 'seqtk subseq';
target_files=(`grep $uo <(echo ${transcripts[@]} | sed 's/ /\n/g') `) ;
echo ${target_files[@]} > "$4"_uotargets.txt ;
unzip -p $2 ${target_files[@]} | awk 'BEGIN{RS=">"} {print $1}' | shuf -n 1000 - | tee "$4"_sampled_seqs.txt | seqtk subseq $3 - | tee "$cat_readfile" > /dev/null ;
samps=$(wc -l "$4"_sampled_seqs.txt | cut -d " " -f 1); if [ $samps -lt 1000 ]
	then printf "$uo\t$samps\n" >> $4_null_samples.txt 
fi ;
#done > $cat_readfile

#fetch matching fast5 files
echo 'fast5 fetcher';

python3 ~/TEST_squiggle/SquiggleKit/fast5_fetcher_multi.py -q $cat_readfile  \
-s $(dirname $3)/sequencing_summary* -o $tmp_dir -m $1;

echo 'tombo preprocess';
tombo preprocess annotate_raw_with_fastqs --overwrite --fast5-basedir $tmp_dir --fastq-filenames $cat_readfile ;

echo 'tombo resquiggle';
cat `find $consensus_ref -name "*${uo};end*.ref.fa" | sort | head -n 1` | tee $4_current_reference.fa >> $4_tombo_reference.fa;
#ref coords are wrong
#grep '>' $4_tombo_reference.fa | cut -d ':' -f 2 | sort | uniq > ref_coords.txt

tombo resquiggle --rna --overwrite --ignore-read-locks $tmp_dir $4_current_reference.fa;
mv $tmp_dir/* $raw_dir;
done
rm -d $tmp_dir
#END LOOP

echo 'tombo detect modifications'
tombo detect_modifications alternative_model --rna --fast5-basedirs $raw_dir --statistics-file-basename $4 --alternate-bases 5mC

echo 'tombo text output'
tombo text_output browser_files --browser-file-basename $4 --file-types dampened_fraction  \
--statistics-filename $4.5mC.tombo.stats --fast5-basedirs $raw_dir
