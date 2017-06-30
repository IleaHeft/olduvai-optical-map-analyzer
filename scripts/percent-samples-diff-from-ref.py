#! /usr/bin/env python

import sys
import pdb
from collections import defaultdict
from collections import Counter


sv_calls = sys.argv[1]

samples_with_diff = defaultdict(list)
num_samples_for_gene = defaultdict(list)

for line in open(sv_calls):
    if line.startswith("#"):
        continue
    else:
        fields = line.strip().split("\t")

        gene = fields[0]
        sample = fields[1]
        diff_from_ref = abs(float(fields[3]))
        
        gene_sample = gene + "_" + sample

        if diff_from_ref >= 2000:
            samples_with_diff[gene].append(sample)

            num_samples_for_gene[gene].append(sample)             
        else:
            num_samples_for_gene[gene].append(sample)             

for key, value in samples_with_diff.items():
    samples = set(value)
    print key,"\t",len(samples)


        
