

### Run _runSV.py_ on each sample

As of 3/26/14, when I ran the SV pipeline on HI1711A and SL1954, I ran it with the "commerical pipeline"
available at the time.  I think all of the 4 1000 Genomes Samples I looked at so far were run with the more 
advanced "in development" pipeline.  Perhaps that could explain why we don't see the supposedly "non-pathogenic"
variation in NBPF1 and NBPF12 in either the schizophrenia or autism samples - but could also just be because
we only have 1 sample for each, so could just be a sample size thing.  Also, my analysis of SL1954 does NOT
shows the 30kb inseretion within NBPF20 like is shown in the image from Pui (even though Pui's was supposedly
also aligned to hg38) - so need to follow up on that.  Need to ask Annie for (1) any additional SV files on 1000
Genomes Samples (and confirm SV pipeleine used), and (2) for their SV call files that they used for the autism and 
schizophrenia samples - and also confirm same pipeline used or not.  Also, need to follow up with Ahmed 
on what exactly is going on/how to interpret the bed file SV output in terms of insertion & delection.  Also,
just need to get IrysView working!

``` bsub < run-sv-on-HI1711A.sh ```

``` bsub < run-sv-on-SL1954.sh ```


### Generate custom track files based on the SV results
#### For a single file

Run the scrip _process-sv-beds-for-custom-track-single-sample.sh_ providing the following arguments:  
1. the SV bed file  
2. the sample ID  
3. the results folder
4. the assembly version (e.g. hg38)  

Code for when I ran this for HI1711A


``` bash scripts/process-sv-beds-for-custom-track-single-sample.sh results/sv-output-HI1711A/merged_smaps/exp_refineFinal1_merged_filter.bed HI1711A results/custom-tracks/ hg38 ```

Code for when I ran this for SL1954

``` bash scripts/process-sv-beds-for-custom-track-single-sample.sh results/sv-output-SL1954/merged_smaps/exp_refineFinal1_merged_filter.bed SL1954 results/custom-tracks/ hg38 ```

#### For a set of files


### Intersect the SV results with a list of genes of interest (e.g. DUF1220 containing) to report out any overlaps.  
Also run the interesect with a certain window of base pairs on either side of the genes of interest

#### For a single sample

Run script _intersect-sv-bed-file-with-coords-single-sample.sh_

Inputs:
1. The SV bed file  
2. The sample id
3. The file giving the coordinates of the genes of interest (already sorted for intersection)  
4. The amount of slop you want to add to either side of each gene  
5. The genome file that specifies the size of chromosomes (necessary so that you don't add slop that runs off the edge of a chromosomes)  
6. The results folder  

``` bash scripts/intersect-sv-bed-file-with-coords-single-sample.sh results/sv-output-SL1954/merged_smaps/exp_refineFinal1_merged_filter.bed SL1954 supporting-reference-data/duf1220-containing-gene-cords-combined-annotations.bed 10000 supporting-reference-data/hg38.genome results/intersections-with-gene-coords/schizophrenia/ ```

``` bash scripts/intersect-sv-bed-file-with-coords-single-sample.sh results/sv-output-HI1711A/merged_smaps/exp_refineFinal1_merged_filter.bed HI1711A supporting-reference-data/duf1220-containing-gene-cords-combined-annotations.bed 10000 supporting-reference-data/hg38.genome results/intersections-with-gene-coords/autism/ ```

#### For multiple samples
