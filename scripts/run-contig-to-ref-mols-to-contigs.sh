#! /usr/bin/env bash

#BSUB -J mols-contigs-ref
#BSUB -o logs/mols-contigs-ref%J.out
#BSUB -e logs/mols-contigs-ref%J.err

source scripts/config.sh

# Provide the "generic" file prefix
generic_prefix=local_contigs_to_ref

region_to_measure=$1
echo $region_to_measure

# Run the contig to ref analysis

contig_to_ref=$output_dir/$num_samples-contigs-to-ref-$region_to_measure-region.txt
echo $contig_to_ref

rm $contig_to_ref

for sample in $(ls $sample_dir | less | sed 's/\///g');
    do
        
        echo "generating duf nick data for: " $sample
        
        duf_nicks=$duf_nick_dir/duf-nicks-$sample.txt
        grep -v "#" $sample_dir/$sample/local_contigs_to_ref_r.cmap | cut -f 5,6 | sed 's/\.[0-9]//g' | awk 'BEGIN{OFS="\t"} {print "chr"wq$1,$2,$2+1}' | sort -k 1,1 -k 2,2n | grep -v "chr0" | bedtools intersect -wa -wb -a stdin -b ~/LabProjects/Irys/annotation-clade-based-numbering-full-domains-2016-11-29.bed > $duf_nicks

        echo "running python script for: " $sample
        python scripts/contigs-to-ref.py $sample_dir/$sample $shift_nicks $region_to_measure $duf_nicks $output_dir $sample $generic_prefix >> $contig_to_ref
    done


# Run the molecule to contig analysis
output_file=$output_dir/$num_samples-mols-to-contigs-$region_to_measure-region.txt
rm $output_file

python scripts/mols-to-contigs.py $sample_dir $shift_nicks $region_to_measure . exp_refineFinal1_contig $contig_to_ref > $output_file
