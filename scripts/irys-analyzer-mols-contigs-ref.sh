#! /usr/bin/env bash
#BSUB -J irys_analyze
#BSUB -o logs/irys_analyze_%J.out
#BSUB -e logs/irys_analyze_%J.err

source scripts/config.sh


rm $output_dir/status.txt
rm $output_dir/ref-distances.tab

touch $output_dir/status.txt
touch $output_dir/ref-distances.tab


#Run the contig to ref and mol to cotig analysis for the HLS region
echo "calculating contig to ref and mol to contig distances for the HLS region"
bash scripts/run-contig-to-ref-mols-to-contigs.sh HLS


#Run the contig to ref and mol to cotig analysis for the HLS region
echo "calculating contig to ref and mol to contig distances for the CON1 region"
bash scripts/run-contig-to-ref-mols-to-contigs.sh CON1

hls_output=$output_dir/$num_samples-mols-to-contigs-HLS-region.txt
con1_output=$output_dir/$num_samples-mols-to-contigs-CON1-region.txt

# run the peak caller

echo "calling HLS region peaks for" $sample
bash $script_dir/peak-caller.sh HLS $hls_output $link_dist $output_dir

echo "calling CON1 region peaks for" $sample
bash $script_dir/peak-caller.sh CON1 $con1_output $link_dist $output_dir


# run the sv caller

echo "calling SVs for the HLS region for" $sample
python $script_dir/sv-caller.py $ref_dist_hls $peak_calls_hls HLS $link_dist $output_dir

echo "calling SVs for the CON1 region for" $sample
python $script_dir/sv-caller.py $ref_dist_con1 $peak_calls_con1 CON1 $link_dist $output_dir

# run the zygosity caller

echo "calling zygosity for the HLS region for" $sample
python $script_dir/call-zygote-status.py $output_dir/sv-calls-HLS-$link_dist.txt $link_dist $output_dir $min_mols_in_cluster


echo "calling zygosity for the CON1 region for" $sample
python $script_dir/call-zygote-status.py $output_dir/sv-calls-CON1-$link_dist.txt $link_dist $output_dir $min_mols_in_cluster


# Calculate the number of structural alleles for each gene in the population analyzed
echo "calculating number of size allels per gene for HLS region for" $sample
bash $script_dir/allele-caller.sh $output_dir/sv-calls-HLS-$link_dist.txt $min_dist $ref_dist_hls > $output_dir/allele-counts-by-gene-hls.txt


echo "calculating number of size allels per gene for CON1 region for" $sample
bash $script_dir/allele-caller.sh $output_dir/sv-calls-CON1-$link_dist.txt $min_dist $ref_dist_con1 > $output_dir/allele-counts-by-gene-con1.txt


# generate files where I have filtered out molecules where the nicks of interest are not adjacent to one another (e.g. there is a nick in the molecule between the aligned CON2 and CON3 nicks)
echo "filtering out molecules without adjacent nicks from HLS results"
awk 'BEGIN{OFS = "\t"} {if ($5 == 1) print $0}' $hls_output > $hls_output_adj_only

echo "filtering out molecules without adjacent nicks from HLS results"
awk 'BEGIN{OFS = "\t"} {if ($5 == 1) print $0}' $con1_output > $con1_output_adj_only


###############################
# execute the peak calling, sv calling, and zygosity calling on only molecules with adjacent nicks

# run the peak caller

echo "calling HLS region peaks for" $sample
bash $script_dir/peak-caller.sh HLS $hls_output_adj_only $link_dist $output_dir_adjonly

echo "calling CON1 region peaks for" $sample
bash $script_dir/peak-caller.sh CON1 $con1_output_adj_only $link_dist $output_dir_adjonly


# run the sv caller

echo "calling SVs for the HLS region for" $sample
python $script_dir/sv-caller.py $ref_dist_hls $peak_calls_hls_adjonly HLS $link_dist $output_dir_adjonly

echo "calling SVs for the CON1 region for" $sample
python $script_dir/sv-caller.py $ref_dist_con1 $peak_calls_con1_adjonly CON1 $link_dist $output_dir_adjonly

# run the zygosity caller

echo "calling zygosity for the HLS region for" $sample
python $script_dir/call-zygote-status.py $output_dir_adjonly/sv-calls-HLS-$link_dist.txt $link_dist $output_dir_adjonly $min_mols_in_cluster


echo "calling zygosity for the CON1 region for" $sample
python $script_dir/call-zygote-status.py $output_dir_adjonly/sv-calls-CON1-$link_dist.txt $link_dist $output_dir_adjonly $min_mols_in_cluster

# Calculate the number of structural alleles for each gene in the population analyzed
echo "calculating number of size allels per gene for HLS region for" $sample
bash $script_dir/allele-caller.sh $output_dir_adjonly/sv-calls-HLS-$link_dist.txt $min_dist $ref_dist_hls > $output_dir_adjonly/allele-counts-by-gene-hls.txt


echo "calculating number of size allels per gene for CON1 region for" $sample
bash $script_dir/allele-caller.sh $output_dir_adjonly/sv-calls-CON1-$link_dist.txt $min_dist $ref_dist_con1 > $output_dir_adjonly/allele-counts-by-gene-con1.txt
