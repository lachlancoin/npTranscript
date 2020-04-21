#!/usr/bin/env python3

import pysam
import numpy
import argparse
import sys

#parse arguments
parser = argparse.ArgumentParser(description = 'Extract alignmened segments of reads from bam file')
required_args = parser.add_argument_group('Positional arguments')
required_args.add_argument('BAM_in', type = str, help = 'Your minibam file containing reads mapping to a certain region')


settings_args = parser.add_argument_group('Settings')
settings_args.add_argument('--jaspar', type = str, help = 'Name of out file to produce count matrix in jaspar .cm format', default = None)
settings_args.add_argument('--fasta', type = str, help = 'Name of out file to produce fasta output without gaps', default = None)
settings_args.add_argument('--coords', metavar=('START','STOP'),required = True, type = int, nargs = 2, help = 'Start and stop co-ordinates for your region of interest')
settings_args.add_argument('--chromosome', type = str, help = 'The chromosome/contig you are using as a reference [Default=None]', default = None)


args = parser.parse_args()

if args.fasta == None and args.jaspar == None:
	sys.exit('Error: Need at least one of --jaspar or --fasta to be selected')

args.coords = [c-1 for c in args.coords] #correct for 0-based indexing
#-------------main------------

#open fasta file for continual writing
if args.fasta:
	f = open(args.fasta, 'w')

sam = pysam.AlignmentFile(args.BAM_in)


motif_seqs = []
motif_length = args.coords[1] - args.coords[0] +1 #window size of all bases in motif
counts = numpy.zeros((4,motif_length))

seq_i = 0

for a in sam.fetch(contig = args.chromosome, start = args.coords[0], stop = args.coords[1]):
	# ref_motif_coord = 64-a.reference_start
	# query_motif_coord = ref_motif_coord + a.query_alignment_start
	# motif_seqs.append(a.query_sequence[query_motif_coord: query_motif_coord + 11])
	seq_i+=1
	pairs = a.get_aligned_pairs()
	pairs_iter = iter(pairs)
	motif_segment = ''
	fasta_segment = ''
	for p in pairs:
		try:
			if p[1] > args.coords[1] +1:
				break
		except TypeError:
			continue
		if p[1] >= args.coords[0] and p[1] <= args.coords[1]:
			if p[0] != None:
				if args.jaspar:
					motif_segment+=a.query_sequence[p[0]]
				if args.fasta:
					fasta_segment+=a.query_sequence[p[0]]
			elif args.jaspar:
				motif_segment+='-'
	if args.jaspar:
		if len(motif_segment) == motif_length:
			print(motif_segment)
			for i, sym in enumerate(motif_segment):
				if sym == 'A':
					counts[0,i] +=1
				if sym == 'C':
					counts[1,i]+=1
				if sym == 'G':
					counts[2,i]+=1
				if sym == 'T':
					counts[3,i]+=1
	if args.fasta:
		#limit to 100 characters per line
		if fasta_segment:
			f.write('>align-region-'+str(seq_i)+'\tcoords:'+'-'.join([str(x) for x in args.coords])+'\t'+a.query_name+'\n')
			for i in range(0, len(fasta_segment), 100):
				try:
					f.write(fasta_segment[i:i+100]+'\n')
				except IndexError:
					f.write(fasta_segment[i::]+'\n')

if args.fasta:
	f.close()

if args.jaspar:
	with open(args.jaspar, 'w') as j:
		for line, sym in zip(counts, ['A| ','C| ','G| ','T| ']):
			count_string = ' '.join([str(i) for i in line])
			j.write(sym+count_string+'\n')
