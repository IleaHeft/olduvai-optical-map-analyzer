contig_to_ref=results/contig-to-ref/3trios-contig-to-ref.txt
data=data/1704_trios_mols_to_contigs_to_ref/mols_to_contigs_to_ref_trios/
output_dir=results
output_file=$output_dir/3trios-mols-to-contigs.txt

python scripts/mols-to-contigs.py $data 0 HLS . exp_refineFinal1_contig $contig_to_ref > $output_file
