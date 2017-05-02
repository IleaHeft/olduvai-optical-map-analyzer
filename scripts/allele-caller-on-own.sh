output_dir=$1
ref_dist_con1=$2
ref_dist_hls=$3
min_dist=$4

bash scripts/allele-caller.sh $output_dir/sv-calls-CON1-2000.txt $min_dist $ref_dist_con1 > $output_dir/allele-counts-by-gene-con1.txt

bash scripts/allele-caller.sh $output_dir/sv-calls-HLS-2000.txt $min_dist $ref_dist_hls > $output_dir/allele-counts-by-gene-hls.txt
