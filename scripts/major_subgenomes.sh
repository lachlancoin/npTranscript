#!/bin/bash -e

#extract lengths of reads with leader
seqtk subseq $READS <(samtools view leader_mapped.bam | cut -f 1 -| sort | uniq) > leader_reads.fq
seqtk comp leader_reads.fq | tee read_comp.txt | cut -f 2 > leader_read_lengths.txt

#R script to find main subgenome mRNAs
./find_subgenomes.R leader_read_lengths.txt > bin_bounds.txt

#get reads for each length bin
#TODO - change temp_file to pipe?

i=1                                                                                                                     
while read line; do bounds=( $line ); file_name="bin${i}_${bounds[0]}_${bounds[1]}.fq" ; \                              
temp_file=$(mktemp); awk -v min=${bounds[0]} -v max=${bounds[1]} '{if ($2 >= min && $2 <= max) print $1}' read_comp.txt > $temp_file; \
seqtk subseq leader_reads.fq $temp_file > $file_name; rm $temp_file; i=$(expr $i + 1); \
done < bin_bounds.txt