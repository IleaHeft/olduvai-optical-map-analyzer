#! /usr/bin/env bash


sample=HG02490
#xmap=data/g4-labeled/xmap-filt1/$sample-filt1.xmap
#xmap_filt_to_display=for-align-vis/$sample-xmap-filt-all-conditions-g4labeled.xmap
#data_file=results/MolRef/HG02490-G4samples-link1000-minmols1-filtmulti1-20170630/con2-con3-dist-all-samples-plus-0.txt



#sample=GM19921
xmap=data/1704_alignments_multimatch_154complete/xmap-filt1/$sample-filt1.xmap
xmap_filt_to_display=for-align-vis/$sample-xmap-filt-all-conditions-standard.xmap
data_file=results/MolRef/154samples-link1000-minmols1-filtmulti1-20170613/con2-con3-dist-all-samples-plus-0.txt


grep $sample $data_file | cut -f 3 > for-align-vis/molecules-to-display-$sample.txt

python scripts/filter-xmap-based-on-mols-that-pass-all-criteria.py for-align-vis/molecules-to-display-$sample.txt $xmap > $xmap_filt_to_display 

