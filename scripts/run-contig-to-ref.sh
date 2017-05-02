#! /usr/bin/env bash

#BSUB -J mols-contigs-ref
#BSUB -o mols-contigs-ref%J.out
#BSUB -e mols-contigs-ref%J.err

# Provide the directory that contains the directory for each sample
sample_dir=data/1704_mols_to_contigs_to_ref_154complete
# Provide the "generic" file prefix
generic_prefix=local_contigs_to_ref

# amount to shift the nicks
shift_nick=0

# region to measure
region_to_measure=HLS

output_dir=results/contig-to-ref
output_file=$output_dir/153samples-contig-to-ref.txt

rm $output_file

for sample in $(ls $sample_dir | less | sed 's/\///g');
    do
        
        echo "generating duf nick data for: " $sample
        
        duf_nicks=duf-nicks-$sample.txt
        grep -v "#" $sample_dir/$sample/local_contigs_to_ref_r.cmap | cut -f 5,6 | sed 's/\.[0-9]//g' | awk 'BEGIN{OFS="\t"} {print "chr"wq$1,$2,$2+1}' | sort -k 1,1 -k 2,2n | grep -v "chr0" | bedtools intersect -wa -wb -a stdin -b ~/LabProjects/Irys/annotation-clade-based-numbering-full-domains-2016-11-29.bed > $duf_nicks

        echo "running python script for: " $sample
        python scripts/contigs-to-ref.py $sample_dir/$sample $shift_nick $region_to_measure $duf_nicks $output_dir $sample $generic_prefix >> $output_file
    done
