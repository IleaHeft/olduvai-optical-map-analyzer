#! /usr/bin/env bash

# this automatically generates todays date
date=$(date +"%Y%m%d")

# the data folder (that can contain different sub-folders of different sample groups)
data_dir=/vol7/home/eskildseni/LabProjects/irys-duf1220/data
# the specific folder that holds all of the samples you want to analyze
sample_dir=$data_dir/g4-labeled


# specify the alignment type (MolRef, MolContig, ContigRef)
alignment_type=MolRef

# number of samples
num_samples=HG02490-G4

# The minimum difference that you want between the maximum confidence score for a molecule and the next highest confidence score - molecules with confidence scores closer than this threshold will be filtered out
conf_spread=1

# The number of base pairs within which to string multiple molecules together in the peak caller step
link_dist=1000

# the file name format for the reference cmap files: everything that follows the sample name
#file_exten_generic=duf1220_mols_v_chr1
#file_exten_rcmap=duf1220_mols_v_chr1_r.cmap
#file_exten_qcmap=duf1220_mols_v_chr1_q.cmap
#file_exten_xmap=duf1220_mols_v_chr1.xmap

file_exten_generic=DufG4_duf1220_new_mols_v_chr1grna_M1
file_exten_rcmap=DufG4_duf1220_new_mols_v_chr1grna_M1_r.cmap
file_exten_qcmap=DufG4_duf1220_new_mols_v_chr1grna_M1_q.cmap
file_exten_xmap=DufG4_duf1220_new_mols_v_chr1grna_M1.xmap


# When generating zygosity calls, specify the minimum number of molecules in each cluster that you want to "count" the SV call for 
## For example, if you only want to consider SV calls supported by 2 or more molecules, set the value to 2 - if you're fine with SV calls only having the support of a single molecule
## then a setting of 1 is fine.  Could also set it higher if you wanted greater levels of support
min_mols_in_cluster=1

results_folder=$alignment_type/$num_samples"samples-link"$link_dist"-minmols"$min_mols_in_cluster"-filtmulti"$conf_spread"-"$date

# DUF annotation file
duf_annotation=~/LabProjects/Irys/annotation-clade-based-numbering-full-domains-2016-11-29.bed

# folder containing the scripts
script_dir=~/LabProjects/irys-duf1220/scripts



# you can modify this to change the CON3 nick to a nick to the right of CON3 by however many nicks you specify - might be useful for seeing if different nicks recover more molecules
shift_nicks=0

# generate directory for the filtered xmaps to go in
filt_xmap_dir=$sample_dir/"xmap-filt"$conf_spread

if [ ! -d "$filt_xmap_dir" ]; then
    mkdir -p $filt_xmap_dir
fi

# the path to the desired output folder
output_dir=/vol7/home/eskildseni/LabProjects/irys-duf1220/results/$results_folder
if [ ! -d "$output_dir" ]; then
    mkdir -p $output_dir
fi

# the path to the desired output folder with a subfolder that contains the results for only the molecules with adjacent nicks
output_dir_adjonly=$output_dir/adj-mols-only

if [ ! -d "$output_dir_adjonly" ]; then
    mkdir -p $output_dir_adjonly
fi

# the names of the output files that contain each relevant molecule and the distance between the relevant nicks on each molecule
hls_output=$output_dir/con2-con3-dist-all-samples-plus-$shift_nicks.txt
con1_output=$output_dir/con1-region-dist-all-samples.txt

hls_output_adj_only=$output_dir_adjonly/con2-con3-dist-all-samples-plus-$shift_nicks.txt
con1_output_adj_only=$output_dir_adjonly/con1-region-dist-all-samples.txt

# the name of the directory and the file that will be created (and overwritten) for each sample that uses bedtools intersect to create a BED file of nick sites within DUF1220 domains
duf_nick_dir=$output_dir/duf-nick-data

if [ ! -d $duf_nick_dir ]; then
    mkdir -p $duf_nick_dir
fi

duf_nicks=$duf_nick_dir/duf-nick-sites.bed


ref_dist_hls=$output_dir/ref-dist-hls.tab
ref_dist_con1=$output_dir/ref-dist-con1-region.tab

peak_calls_hls=$output_dir/HLS-peak-calls-$link_dist.txt
peak_calls_con1=$output_dir/CON1-peak-calls-$link_dist.txt

peak_calls_hls_adjonly=$output_dir_adjonly/HLS-peak-calls-$link_dist.txt
peak_calls_con1_adjonly=$output_dir_adjonly/CON1-peak-calls-$link_dist.txt
