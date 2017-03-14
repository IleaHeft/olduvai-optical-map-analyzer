#! /usr/bin/env bash

date=$(date +"%Y%m%d")
#today=$(date +"%Y-%m-%d")
#today=2017-03-02
today=trios-T8A5-link2000-min2mols-filtmulti-$date

# the folder that holds all of the samples
sample_dir=/vol7/home/eskildseni/LabProjects/irys-duf1220/data/trios-T8A5-multimatch

# folder containing the scripts
script_dir=~/LabProjects/irys-duf1220/scripts


# the file name format for the reference cmap files: everything that follows the sample name
file_exten_generic=duf1220-mols-v-chr1
file_exten_rcmap=duf1220-mols-v-chr1-r.cmap
file_exten_qcmap=duf1220-mols-v-chr1-q.cmap
file_exten_xmap=duf1220-mols-v-chr1.xmap

# The minimum difference that you want between the maximum confidence score for a molecule and the next highest confidence score - molecules with confidence scores closer than this threshold will be filtered out
conf_spread=1

# The number of base pairs within which to string multiple molecules together in the peak caller step
link_dist=2000

# When generating zygosity calls, specify the minimum number of molecules in each cluster that you want to "count" the SV call for 
## For example, if you only want to consider SV calls supported by 2 or more molecules, set the value to 2 - if you're fine with SV calls only having the support of a single molecule
## then a setting of 1 is fine.  Could also set it higher if you wanted greater levels of support
min_mols_in_cluster=2


# you can modify this to change the CON3 nick to a nick to the right of CON3 by however many nicks you specify - might be useful for seeing if different nicks recover more molecules
shift_nicks=0


# the path to the desired output folder
output_dir=/vol7/home/eskildseni/LabProjects/irys-duf1220/results/$today
mkdir $output_dir

# the path to the desired output folder with a subfolder that contains the results for only the molecules with adjacent nicks
output_dir_adjonly=$output_dir/adj-mols-only
mkdir $output_dir_adjonly

# the names of the output files that contain each relevant molecule and the distance between the relevant nicks on each molecule
hls_output=$output_dir/con2-con3-dist-all-samples-plus-$shift_nicks.txt
con1_output=$output_dir/con1-region-dist-all-samples.txt

hls_output_adj_only=$output_dir_adjonly/con2-con3-dist-all-samples-plus-$shift_nicks.txt
con1_output_adj_only=$output_dir_adjonly/con1-region-dist-all-samples.txt

# the name of the file that will be created (and overwritten) for each sample that uses bedtools intersect to create a BED file of nick sites within DUF1220 domains
duf_nicks=$output_dir/duf-nick-sites.bed


ref_dist_hls=$output_dir/ref-dist-hls.tab
ref_dist_con1=$output_dir/ref-dist-con1-region.tab

peak_calls_hls=$output_dir/HLS-peak-calls-$link_dist.txt
peak_calls_con1=$output_dir/CON1-peak-calls-$link_dist.txt

peak_calls_hls_adjonly=$output_dir_adjonly/HLS-peak-calls-$link_dist.txt
peak_calls_con1_adjonly=$output_dir_adjonly/CON1-peak-calls-$link_dist.txt
