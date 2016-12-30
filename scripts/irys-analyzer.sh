#! /usr/bin/env bash
#BSUB -J nick_dist
#BSUB -o nick_dist_%J.out
#BSUB -e nick_dist_%J.err
#BSUB -n 12

source scripts/config.sh

#sample_dir=$1
#shift_nicks=$2





rm $output_dir/status.txt
rm $output_dir/ref-distances.tab
rm $hls_output
rm $con1_output

touch $output_dir/status.txt
touch $output_dir/ref-distances.tab
touch $hls_output
touch $con1_output



for sample in $(ls $sample_dir);
    do
        echo $sample
        align_mol_dir=$sample_dir$sample/output/contigs/alignmolvref/merge/
        
        contig1_rcmap=$align_mol_dir/alignmolvref_contig1_r.cmap

        echo $align_mol_dir
    
        grep -v "#" $contig1_rcmap | cut -f 5,6 | sed 's/\.[0-9]//g' | awk 'BEGIN{OFS="\t"} {print "chr"wq$1,$2,$2+1}' | sort -k 1,1 -k 2,2n | grep -v "chr0" |\
        bedtools intersect -wa -wb -a stdin -b ~/LabProjects/Irys/annotation-clade-based-numbering-full-domains-2016-11-29.bed > $duf_nicks
        
        python $script_dir/nick-distance-calc.py $align_mol_dir $shift_nicks HLS $duf_nicks $output_dir >> $hls_output

        python $script_dir/nick-distance-calc.py $align_mol_dir $shift_nicks CON1 $duf_nicks $output_dir >> $con1_output
    done


# run the peak caller
bash $script_dir/peak-caller.sh HLS $hls_output $link_dist

bash $script_dir/peak-caller.sh CON1 $con1_output $link_dist


# run the sv caller
#ref_dist_hls=ref-dist-hls.tab
#ref_dist_con1=ref-dist-con1-region.tab
#
#peak_calls_hls=HLS-peak-calls-$link_dist.txt
#peak_calls_con1=CON1-peak-calls-$link_dist.txt


python $script_dir/sv-caller.py $ref_dist_hls $peak_calls_hls HLS 1000 $output_dir

python $script_dir/sv-caller.py $ref_dist_con1 $peak_calls_con1 CON1 1000 $output_dir

