#! /usr/bin/env bash

sample=HG02490
# you can only sort based on the SV size of one gene at a time
gene=NBPF1

# standard or g4 labeled
label=standard

output_dir=for-align-vis

ordered_xmap=$output_dir/$sample-ordered-$label.xmap

xmap_filt_to_display=$output_dir/$sample-xmap-filt-all-conditions-standard.xmap

data_file=results/MolRef/154samples-link2000-minmols1-filtmulti1-20170611/con2-con3-dist-all-samples-plus-0.txt

order_data=$output_dir/$sample-ordered-mols-$gene.txt


grep $sample $data_file | grep -w $gene | sort -k 9,9n > $order_data

python scripts/order-xmap-for-vis.py $xmap_filt_to_display $order_data > $ordered_xmap

