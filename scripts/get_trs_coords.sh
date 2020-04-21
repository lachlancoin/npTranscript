#!/bin/bash

#MEME SUITE location
MS='/home/lcoin/meme-5.1.1/'
VIRUS=''

#check meme suite is accessible
command -v fimo >/dev/null 2>&1 || { command -v $MS/src/fimo >/dev/null 2>&1; } || { echo >&2 "FIMO not found.  Aborting."; exit 1; }
command -v jaspar2meme >/dev/null 2>&1 || { command -v $MS/scripts/jaspar2meme >/dev/null 2>&1; } || { echo >&2 "jaspar2meme not found.  Aborting."; exit 1; }

#Find TRS. Shift this to other script
#TRS-L coords. Exact string match. Ouputs co-ord file "TRS_short.coords.csv"
python3 trs-coords.py --virus_genome $VIRUS

#coords of first match (ie. TRS-L)
declare -a coord_array=(`sed '2q;d' TRS_short.coords.csv | sed -r 's/[\,|A-Z]/ /g'`)
(( coord_array[0]-=5 )) #extend out to include preceeding conserved hexamer


#build motif. Modify so coords are extracted automatically
python3 motif_from_bam.py --jaspar TRS_long.cm --coords ${coord_array[@]} sorted.host_mapped.bam


$MS/scripts/jaspar2meme -strands 1 -cm .
cat $VIRUS | $MS/src/fimo --norc --oc TRS_long_fimo --thresh 1e-3 trs_only.meme -
