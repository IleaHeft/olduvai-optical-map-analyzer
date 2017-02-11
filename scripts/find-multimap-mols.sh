#! /usr/bin/env bash
#BSUB -J findqmols
#BSUB -o logs/findqmols_%J.out
#BSUB -e logs/findqmols_%J.err

source scripts/config-questionable-mols.sh

duf_annotation=/vol5/home/eskildseni/LabProjects/Irys/annotation-clade-based-numbering-full-domains-2016-11-29.bed 
shift_nicks=0
conf_spread=2
mols_to_strip_hls=$output_dir/mols-to-strip-hls-$conf_spread.txt
mols_to_strip_con1=$output_dir/mols-to-strip-con1-$conf_spread.txt

for sample in $(ls $sample_dir | less | cut -f 1 -d "-" | sort | uniq );
    do
        echo $sample
        align_mol_dir=$sample_dir
        
        contig1_rcmap=$align_mol_dir/$sample-duf1220-mols-v-chr1-r.cmap
        echo $contig1_rcmap
        
    
        grep -v "#" $contig1_rcmap | cut -f 5,6 | sed 's/\.[0-9]//g' | awk 'BEGIN{OFS="\t"} {print "chr"wq$1,$2,$2+1}' | sort -k 1,1 -k 2,2n | grep -v "chr0" |\
        bedtools intersect -wa -wb -a stdin -b $duf_annotation > $duf_nicks
        
        python $script_dir/nick-distance-calc-multi2.py $align_mol_dir $shift_nicks HLS $duf_nicks $output_dir $sample $conf_spread >> $mols_to_strip_hls

        python $script_dir/nick-distance-calc-multi2.py $align_mol_dir $shift_nicks CON1 $duf_nicks $output_dir $sample $conf_spread >> $mols_to_strip_con1
    done



