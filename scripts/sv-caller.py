#! /usr/bin/env python

import sys
import pdb
from collections import defaultdict

ref_distances=sys.argv[1]
peak_calls=sys.argv[2]
region = sys.argv[3]
link_dist = sys.argv[4]

# set up output file
out = open("sv-calls-" + region + "-" + link_dist + ".txt", mode = 'w')

# set up output file for unused molecules that don't pass the current filter of needing at least two molecules in a bin
unused_mols = open("unused-mols-sv-caller.txt",mode = 'w')

list_of_samples=[]
list_of_genes=[]
# Generate a dictionary of genes and the distance in the reference between nicks

#initialize dictionary
# keys: gene
# value: distance between nicks in the reference

gene_ref_dist=defaultdict(int)

for line in open(ref_distances):
    fields = line.strip().split("\t")
    sample = fields[0]
    gene = fields[1]
    ref_dist = float(fields[4])
    
    sample_gene = sample + "_" + gene
    gene_ref_dist[sample_gene] = ref_dist


# Do the peak calling

#initialize dictionary 
# keys: sample/gene name combo
# values: diff from reference for each cluster

sample_gene_diffs=defaultdict(list)
gene_nick_type = defaultdict(str)

for line in open(peak_calls):
    fields = line.strip().split("\t")
    sample = fields[0]
    gene = fields[1]
    min_dist = int(fields[2])
    max_dist = int(fields[3])
    mean_dist = float(fields[4])
    median_dist = float(fields[5])
    count = int(fields[6])
    median_conf = fields[7]
    conf_values = fields[8]
    dist_bw_nick_nums=fields[9]
    align_qual = fields[10]
    num_left_nicks = fields[11]
    num_right_nicks = fields[12]
    distances = fields[13]
    mol_ids = fields[14]
    nick_type = fields[15]

    gene_nick_type[gene] = nick_type

    list_of_samples.append(sample)
    list_of_genes.append(gene)

    dist_range = str(abs(max_dist - min_dist))
    sample_gene = sample + "_" + gene

    ref_dist = gene_ref_dist.get(sample_gene)

    diff = ref_dist - median_dist
    
    sample_gene = sample + "_" + gene
    

    if count >= 2:
        

        sample_gene_diffs[sample_gene].append((diff,count,str(min_dist),str(max_dist),dist_range,median_conf))

    # if there is a molecule not binned with other molecules, consider it if it meets certain criteria
    elif count == 1:
        
        # set the conditions under which you'd accept a lone molecule
        # align_qual == "GG" means that the molecule has at least one nick to the left and one nick to the right of the nicks of interest and those flanking nicks are both aligned "in order" - no nicks are skipped on either the reference or the query
        if "U" not in align_qual:
            
            sample_gene_diffs[sample_gene].append((diff,count,str(min_dist),str(max_dist),dist_range,median_conf))
        else:
            
            print >> unused_mols, line.strip()

    else:
        print >> unused_mols, line.strip()


# For each cluster of molecules (or for single molecules that look well aligned), calculate the difference between the distance for that molecule and the reference distance put it in a category based on the difference
for sample_gene, diffs_counts in sample_gene_diffs.items():
    
    for item in diffs_counts:
        
        diff = item[0]
        count = item[1]
        min_dist = item[2]
        max_dist = item[3]
        dist_range = item[4]
        median_conf = item[5]

        if abs(diff) < 2000:
            diff_status = "ref"
        elif diff > 0: 
            
            diff_status = "del"
        elif diff < 0:

            diff_status = "ins"
        else:
            continue

        sample = sample_gene.split("_")[0]
        gene = sample_gene.split("_")[1]

        to_print = [gene,sample,diff_status,str(abs(diff)),str(count),min_dist,max_dist,dist_range,gene_nick_type.get(gene),median_conf]
        print >> out, "\t".join(to_print)

for sample in set(list_of_samples):
    for gene in set(list_of_genes):
        sample_gene = sample + "_" + gene

        if sample_gene not in sample_gene_diffs.keys():
            to_print = [gene, sample, "nd","0","0","0","0","0",gene_nick_type.get(gene),median_conf]
            print >> out, "\t".join(to_print)

