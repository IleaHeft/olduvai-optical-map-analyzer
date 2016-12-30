#! /usr/bin/env bash


today=$(date +"%Y-%m-%d")

# folder containing the scripts
script_dir=~/LabProjects/irys-duf1220/scripts

# the folder that holds all of the samples
sample_dir=/vol5/home/eskildseni/LabProjects/Irys/data/1000GenomesIrys/full-output-folders/

# The number of base pairs within which to string multiple molecules together in the peak caller step
link_dist=1000

# you can modify this to change the CON3 nick to a nick to the right of CON3 by however many nicks you specify - might be useful for seeing if different nicks recover more molecules
shift_nicks=0


# the path to the desired output folder
output_dir=/vol5/home/eskildseni/LabProjects/irys-duf1220/results/$today
mkdir $output_dir

# the names of the output files that contain each relevant molecule and the distance between the relevant nicks on each molecule
hls_output=$output_dir/con2-con3-dist-all-samples-plus-$shift_nicks.tab
con1_output=$output_dir/con1-region-dist-all-samples.tab



# the name of the file that will be created (and overwritten) for each sample that uses bedtools intersect to create a BED file of nick sites within DUF1220 domains
duf_nicks=$output_dir/duf-nick-sites.bed


ref_dist_hls=$output_dir/ref-dist-hls.tab
ref_dist_con1=$output_dir/ref-dist-con1-region.tab

peak_calls_hls=$output_dir/HLS-peak-calls-$link_dist.txt
peak_calls_con1=$output_dir/CON1-peak-calls-$link_dist.txt
