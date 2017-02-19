#! /usr/bin/env bash

#gene=$1
#sample=$2

source scripts/config-yulia-alignments.sh

# region being analyzed
region=$1

# this is the output file from the first nick-distance-calc python script
dist_per_mol=$2

# the link distance is the "bin" size -- how close two molecules have to be together for them to be merged together (e.g 1000 for 1000bp)
link_distance=$3

# this is the output location of the binned/merged data
output_file=$output_dir/$region-peak-calls-$link_distance.txt
#output_all_molecules=merged-dist-$link_distance-all-mols.txt

sample_list=$output_dir/samples.txt
gene_list=$output_dir/genes.txt

# remove any existing version of these files so that you generate a new one and don't just add lines to an old one
rm $output_file
rm $sample_list
rm $gene_list

# generate file with sample names
cut -f 1 $dist_per_mol | sort | uniq > $sample_list

# generate file with gene names
cut -f 2 $dist_per_mol | sort | uniq > $gene_list


# This step merges molecules at similiar positions into "bins" and prints useful information about the molecules in each bin

# Relevant info for understanding bedtools merge operations, is knowing which column of the input file holds what information
# Original input file: $dist_per_mol
# 1. sample
# 2. gene
# 3. molecule ID
# 4. confidence score
# 5. distance between query nick numbers on molecules (generally expect this to be 1 - but a number different than 1 isn't necessarily bad)
# 6. "P" or "G" indicating good or poor alignment at nick sites adjacent to sites of interest
# 7. number of nicks to the left of the left-most nick of interest
# 8. number of nicks to the right of the right-most nick of interest
# 9. the distance between the two nicks of interest

# After the awk line the file looks like a BED file in order to use bedtools merge
# 1. "chr1" - only purpose of this is make bedtools merge work
# 2. The distance between the two nicks of interest
# 3. The distance between the two nicks of interest plus the "link distance" that you input 
# 4. sample
# 5. gene
# 6. molecule ID
# 7. confidence score
# 8. distance between query nick numbers on molecules (generally expect this to be 1 - but a number different than 1 isn't necessarily bad)
# 9. "P" or "G" indicating alignment quality (see explanation above)
# 10. number of nicks to the left of the left-most nick of interest
# 11. number of nicks to the right of the right-most nick of interest
# 12. "std" or "alt" stating whether the coordinates used for the calculation where the standard CON2 and CON3 nicks or whether alternate coordinates were used (because the gene didn't have nicks in both CON2 and CON3)

# The grep -v "PP", etc removes molecules with poor alignments matching those alignment quality markers

# The bedtools merge step returns the following column structure:
# 1. chr1 (only here because of the BED file structure)
# 2. the smallest distance between nicks in the bin
# 3. this is an irrevalent distance - only necessary because of the bedtools merge set-up, it is the largest distance in the bin plus the link distance
# 4. sample name(s) for all molecules in the bin -- there should only be one sample, if there is more than one there is a problem in the code
# 5. gene name for all the molecules in the bin -- there should only be one gene name, if there is more than one there is a problem in the code
# 6. count of the number of molecules in the bin
# 7. list of all confidence scores in the bin
# 8. list of all distances between query nick numbers in the bin
# 9. list of all alignment quality ratings (P or G) in the bin
# 10. list of number of nicks to the left for every molecule in the bin
# 11. list of number of nicks to the right for every molecule in the bin
# 12: list of distances of all molecules in the bin
# 13. list of all the molecules IDs in the bin
# 14. Whether the gene used standard "std" or alternate "alt" nicks

# The final file has columns as specified immediately above, except the first three columns have been removed


for i in $(less $sample_list);
    do
        sample=$i
       # echo $sample

        for i in $(less $gene_list);
            do
                gene=$i
               # echo $gene
            
                awk -v link="$link_distance" 'BEGIN{OFS="\t"} {print "chr1",$9,$9+link,$1,$2,$3,$4,$5,$6,$7,$8,$10}' $dist_per_mol |\
                sed 's/\.[0-9]*//g' |\
                sort -k 1,1 -k 2,2n |\
                grep -w "$gene" |\
                grep $sample |\
                grep -v -P "\-[GU]*P+"\|"\t[GU]*P+" |\
                bedtools merge -i stdin -c 4,5,2,2,2,2,6,7,7,8,9,10,11,2,6,12 -o distinct,distinct,min,max,mean,median,count,median,collapse,collapse,collapse,collapse,collapse,collapse,collapse,distinct |\
                cut -f 4- >> $output_file
            
            
            done

    done




# Make a file of unused molecules
filt_out=$output_dir/unused-mols-bin-step.txt
rm $filt_out

grep -P "\-[GU]*P+"\|"\t[GU]*P+" $dist_per_mol > $filt_out
