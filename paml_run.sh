#!/bin/bash

# first argument should be the protein sequence file.
# second argument should be the cds nucleic acid sequence file.

# alignment of the proteic sequences
clustalw2 -quiet -align -infile=$1 -outfile=prot.ali.aln

# cds sequences aligned by corresponding
./pal2nal.pl prot.ali.aln $2 -output paml > cds.ali.phy

# modification of the control file
awk -v file=cds.ali.phy '{gsub("XXXXX",file); print $0}' yn00.ctl_master > yn00.ctl

#final call, running yn00 program from the paml tool
yn00 > res

#awk pour parser le fichier resultat et choper les KA et KS
