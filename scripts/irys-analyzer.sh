#! /usr/bin/env bash
#BSUB -J irys_analyze
#BSUB -o logs/irys_analyze_%J.out
#BSUB -e logs/irys_analyze_%J.err

source scripts/config.sh


rm $output_dir/status.txt
rm $output_dir/ref-distances.tab
rm $hls_output
rm $con1_output

touch $output_dir/status.txt
touch $output_dir/ref-distances.tab
touch $hls_output
touch $con1_output


echo "sample directory is:" $sample_dir

# Decide which nick distance calculator will be used based on whether you are looking at mol to reference data or mol to contig or contig to ref
if [ $alignment_type == "MolRef" ]; then
    
    echo "alignment_type == MolRef" $alignment_type

    for sample in $(ls $sample_dir | less | cut -f 1 -d "_" | sort | uniq );
        do
            echo $sample
            align_mol_dir=$sample_dir
            
            rcmap=$align_mol_dir/$sample"_"$file_exten_rcmap
            echo "reference cmap file being used is:" $rcmap
            
            echo "generating a bed file of nick sites within DUF1220 domains (start of short exon to end of long exon) for" $sample

            grep -v "#" $rcmap | cut -f 5,6 | sed 's/\.[0-9]//g' | awk 'BEGIN{OFS="\t"} {print "chr"wq$1,$2,$2+1}' | sort -k 1,1 -k 2,2n | grep -v "chr0" |\
            bedtools intersect -wa -wb -a stdin -b $duf_annotation > $duf_nicks
            
            
            # Generate filterted xmap file that removes molecules with secondary mappings of similiar confidence to max confidence and reports only highest confidence alignment for other molecules
            echo "generating a xmap file with multi-match molecules removed for" $sample
            xmap=$sample_dir/$sample"_"$file_exten_xmap
            echo "xmap file being used is:" $xmap
            python $script_dir/generate-filtered-xmap.py $sample_dir $sample $xmap $conf_spread
            
            
            echo "calculating the distance between CON2 and CON3 nicks for" $sample
            python $script_dir/nick-distance-calc.py $align_mol_dir $shift_nicks HLS $duf_nicks $output_dir $sample $file_exten_generic >> $hls_output
            
            echo "calculating the distance between the CON1 nick and the next closest nick upstream for" $sample
            python $script_dir/nick-distance-calc.py $align_mol_dir $shift_nicks CON1 $duf_nicks $output_dir $sample $file_exten_generic >> $con1_output
        done

else

    echo "alignment_type != MolRef, alignment_type is: " $alignment_type
    
    #Run the contig to ref and mol to cotig analysis for the HLS region
    echo "calculating contig to ref and mol to contig distances for the HLS region"
    bash scripts/run-contig-to-ref-mols-to-contigs.sh HLS
    
    #Run the contig to ref and mol to cotig analysis for the HLS region
    echo "calculating contig to ref and mol to contig distances for the CON1 region"
    bash scripts/run-contig-to-ref-mols-to-contigs.sh CON1
    
    if [ $alignment_type == "MolContig" ]; then
        
        echo "using distance files for mols to contig" 
        hls_output=$output_dir/$num_samples-mols-to-contigs-HLS-region.txt
        con1_output=$output_dir/$num_samples-mols-to-contigs-CON1-region.txt

    else
        echo "using distance files for contigs to ref" 
        hls_output=$output_dir/$num_samples-contigs-to-ref-HLS-region.txt
        con1_output=$output_dir/$num_samples-contigs-to-ref-CON1-region.txt
        
        # There is only expected to be one contig for a given region, so makes sense to have the min_mol threshold set to one
        min_mols_in_cluster=1
    fi

fi


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


### make version of the sv all output where homozygous calls are duplicated so that each "allele" is represented by a line
echo "making sv call file with duplicated homozyougs calls for HLS region"
python $script_dir/duplicate-sv-call-if-homozygous.py  $output_dir/sv-calls-HLS-2000.txt > $output_dir/sv-calls-HLS-2000-1lineperallele.txt

echo "making sv call file with duplicated homozyougs calls for CON1 region"
python $script_dir/duplicate-sv-call-if-homozygous.py  $output_dir/sv-calls-CON1-2000.txt > $output_dir/sv-calls-CON1-2000-1lineperallele.txt

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


### make version of the sv all output where homozygous calls are duplicated so that each "allele" is represented by a line
echo "making sv call file with duplicated homozyougs calls for HLS region"
python $script_dir/duplicate-sv-call-if-homozygous.py  $output_dir_adjonly/sv-calls-HLS-2000.txt > $output_dir_adjonly/sv-calls-HLS-2000-1lineperallele.txt

echo "making sv call file with duplicated homozyougs calls for CON1 region"
python $script_dir/duplicate-sv-call-if-homozygous.py  $output_dir_adjonly/sv-calls-CON1-2000.txt > $output_dir_adjonly/sv-calls-CON1-2000-1lineperallele.txt


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
