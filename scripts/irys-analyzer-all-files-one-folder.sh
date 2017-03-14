#! /usr/bin/env bash
#BSUB -J irys_analyze
#BSUB -o logs/irys_analyze_%J.out
#BSUB -e logs/irys_analyze_%J.err

source scripts/config-all-files-one-folder.sh


rm $output_dir/status.txt
rm $output_dir/ref-distances.tab
rm $hls_output
rm $con1_output

touch $output_dir/status.txt
touch $output_dir/ref-distances.tab
touch $hls_output
touch $con1_output



for sample in $(ls $sample_dir | less | cut -f 1 -d "-" | sort | uniq );
    do
        echo $sample
        align_mol_dir=$sample_dir
        
        contig1_rcmap=$align_mol_dir/$sample-$file_exten_rcmap
        echo $contig1_rcmap
        
        echo "generating a bed file of nick sites within DUF1220 domains (start of short exon to end of long exon) for" $sample

        grep -v "#" $contig1_rcmap | cut -f 5,6 | sed 's/\.[0-9]//g' | awk 'BEGIN{OFS="\t"} {print "chr"wq$1,$2,$2+1}' | sort -k 1,1 -k 2,2n | grep -v "chr0" |\
        bedtools intersect -wa -wb -a stdin -b ~/LabProjects/Irys/annotation-clade-based-numbering-full-domains-2016-11-29.bed > $duf_nicks
        
        
        # Generate filterted xmap file that removes molecules with secondary mappings of similiar confidence to max confidence and reports only highest confidence alignment for other molecules
        echo "generating a xmap file with multi-match molecules removed for" $sample
        xmap=$sample_dir/$sample-$file_exten_xmap
        echo $xmap
        python $script_dir/generate-filtered-xmap.py $sample_dir $sample $xmap $conf_spread
        
        
        echo "calculating the distance between CON2 and CON3 nicks for" $sample
        python $script_dir/nick-distance-calc.py $align_mol_dir $shift_nicks HLS $duf_nicks $output_dir $sample $file_exten_generic >> $hls_output
        
        echo "calculating the distance between the CON1 nick and the next closest nick upstream for" $sample
        python $script_dir/nick-distance-calc.py $align_mol_dir $shift_nicks CON1 $duf_nicks $output_dir $sample $file_exten_generic >> $con1_output
    done



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
