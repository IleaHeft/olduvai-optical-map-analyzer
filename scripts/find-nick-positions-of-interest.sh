#! /usr/bin/env bash

# path to the alignmolvref_contig1_r.cmap
contig1_rcmap=$1

grep -v "#" $contig1_rcmap | cut -f 5,6 | sed 's/\.[0-9]//g' | awk 'BEGIN{OFS="\t"} {print "chr"wq$1,$2,$2+1}' | sort -k 1,1 -k 2,2n | grep -v "chr0" |\
bedtools intersect -wa -wb -a stdin -b ~/LabProjects/Irys/annotation-clade-based-numbering-full-domains-2016-11-29.bed > duf-nick-sites.bed
