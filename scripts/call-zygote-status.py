#! /usr/bin/env python

import sys
import pdb
from collections import Counter
from collections import defaultdict

# provide file of SV calls
sv_calls = sys.argv[1]
link_dist = sys.argv[2]
output_dir = sys.argv[3]
min_num_mols = int(sys.argv[4])

# set up output folder
if "HLS" in sv_calls:
    region = "HLS"
elif "CON1" in sv_calls:
    region = "CON1"
else:
    print "ERROR: UNABLE TO DETERMINE REGION"

out = open(output_dir + "/zygosity-calls" + "-" + region + "-" + link_dist + ".txt",mode = 'w')


# set up counter and dictionary
gene_sample_counter = Counter()
gene_sample_calls=defaultdict(list)

# print colnames
colnames = ["#sample","gene","zygosity","num.calls","size.calls"]
print >> out, "\t".join(colnames)

for line in open(sv_calls):
    if "#" in line:
        continue
    else:
        fields = line.strip().split("\t")
        gene = fields[0]
        sample = fields[1]
        type_call = fields[2]
        size = fields[3]
        num_mols = int(fields[4])
        
        if "nd" in type_call:

            to_print = [sample, gene, "nd","0","nd"]
            print >> out, "\t".join(to_print)        
        
        
        # only looks at SV calls supported by more than on molecule
        elif num_mols >= min_num_mols:
            gene_sample = gene + "_" + sample
            gene_sample_counter[gene_sample] += 1
            gene_sample_calls[gene_sample].append(size)

        else:
            continue
for gene_sample, size_calls in gene_sample_calls.items():
    if len(size_calls) == 1:
        zygosity = "homozygous" + "\t" + str(len(size_calls)) + "\t" + ", ".join(size_calls)
    else:
        zygosity = "heterozygous" + "\t" + str(len(size_calls)) + "\t" + ", ".join(size_calls)

    gene = gene_sample.split("_")[0]
    sample = gene_sample.split("_")[1]
    to_print = [sample,gene,zygosity] 
    print >> out,"\t".join(to_print)
