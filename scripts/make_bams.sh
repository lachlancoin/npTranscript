#!/bin/bash -e
READS=''

#map to whole genome
minimap2 -ax splice -un ../data/SARS-Cov2/wuhan_coronavirus_australia.fasta $READS | samtools view -hF4 -b > whole_genome_mapped.bam

#map to leaderr
minimap2 -ax map-ont -un ../data/SARS-Cov2/leader.fasta $READS | satmools view -hF4 -b > leader_mapped.bam