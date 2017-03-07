#! /usr/bin/env python

import sys
import pdb
from collections import Counter
from collections import defaultdict

# provide file of SV calls
sv_calls = sys.argv[1]

# set up counter and dictionary
gene_sample_counter = Counter()
gene_sample_calls=defaultdict(list)

# print colnames
colnames = ["sample","gene","zygosity","num.calls","size.calls"]
print "\t".join(colnames)

for line in open(sv_calls):
    fields = line.strip().split("\t")
    gene = fields[0]
    sample = fields[1]
    type_call = fields[2]
    size = fields[3]
    num_mols = fields[4]
    
    if "nd" in type_call:

        to_print = [sample, gene, "nd","nd","nd"]
        print "\t".join(to_print)        
        
    # only looks at SV calls supported by only on molecule
    elif num_mols != "1":
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
    print "\t".join(to_print)
