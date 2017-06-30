#! /usr/bin/env bash

input_file=$1
#gene=$2
min_dist=$2

gene_list=genes.txt
ref_dist_list=$3


# print the header line to the file
echo "#gene	region.size.in.ref	num.calls.of.this.size	median.size.in.cluster	min.size.in.cluster	max.size.in.cluster	range.of.diff.from.ref.in.cluster"



# make a list of genes to loop over
grep -v "#" $input_file | cut -f 1 | sort | uniq > genes.txt



for gene in $(less $gene_list);
    do
        # Calculate alleles with amplifications - or just a size change relative to the reference that is greater than zero
        
        ref_dist=$(grep -w $gene $ref_dist_list  | cut -f 5 | sort | uniq | sed 's/\.[0-9]*//g')


        grep -v "#\|nd" $input_file |\
        grep -w $gene |\
        awk -v ref_dist="$ref_dist" -v min_dist="$min_dist" 'BEGIN{OFS="\t"} {if ($5 != 1) print "chr1",ref_dist+$4,ref_dist+$4+min_dist,$1,$4,ref_dist+$4}' |\
        sed 's/\.[0-9]*//g' |\
        sort -k 1,1 -k 2,2n |\
        bedtools merge -i stdin -c 4,4,6,6,6,5 -o distinct,count,median,min,max,collapse |\
        awk -v ref_dist="$ref_dist" 'BEGIN{OFS="\t"} {print $4,ref_dist,$5,$6,$7,$8,$8-$7,$9}'

    done

# Remove the genes list file
rm genes.txt
