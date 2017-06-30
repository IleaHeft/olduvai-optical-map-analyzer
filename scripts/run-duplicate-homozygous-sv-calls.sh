#! /usr/bin/env bash

source scripts/config.sh


### make version of the sv all output where homozygous calls are duplicated so that each "allele" is represented by a line
echo "making sv call file with duplicated homozyougs calls for HLS region"
python $script_dir/duplicate-sv-call-if-homozygous.py  $output_dir/sv-calls-HLS-$link_dist.txt > $output_dir/sv-calls-HLS-$link_dist-1lineperallele.txt

echo "making sv call file with duplicated homozyougs calls for CON1 region"
python $script_dir/duplicate-sv-call-if-homozygous.py  $output_dir/sv-calls-CON1-$link_dist.txt > $output_dir/sv-calls-CON1-$link_dist-1lineperallele.txt



# Run for the adjacent mols only folder
### make version of the sv all output where homozygous calls are duplicated so that each "allele" is represented by a line
echo "making sv call file with duplicated homozyougs calls for HLS region"
python $script_dir/duplicate-sv-call-if-homozygous.py  $output_dir_adjonly/sv-calls-HLS-$link_dist.txt > $output_dir_adjonly/sv-calls-HLS-$link_dist-1lineperallele.txt

echo "making sv call file with duplicated homozyougs calls for CON1 region"
python $script_dir/duplicate-sv-call-if-homozygous.py  $output_dir_adjonly/sv-calls-CON1-$link_dist.txt > $output_dir_adjonly/sv-calls-CON1-$link_dist-1lineperallele.txt


