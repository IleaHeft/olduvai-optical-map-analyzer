#! /usr/bin/env bash

# the folder that holds all of the samples
sample_dir=full-output-folders/

# you can modify this to change the CON3 nick to a nick to the right of CON3 by however many nicks you specify - might be useful for seeing if different nicks recover more molecules
shift_nicks=0

# the name of the file that will be created (and overwritten) for each sample that uses bedtools intersect to create a BED file of nick sites within DUF1220 domains
duf_nicks=duf-nick-sites.bed

# the names of the output files that contain each relevant molecule and the distance between the relevant nicks on each molecule
hls_output=con2-con3-dist-all-samples-plus-$shift_nicks.tab
con1_output=con1-region-dist-all-samples.tab

# The number of base pairs within which to string multiple molecules together in the peak caller step
link_dist=1000
