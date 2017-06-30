#! /usr/bin/env python

import sys
import pdb
from collections import defaultdict
from collections import Counter


sv_calls = sys.argv[1]

gene_counts = Counter()
homozygous = []

for line in open(sv_calls):
    if line.startswith("#"):
        continue
    else:
        fields = line.strip().split("\t")
        
        gene = fields[0]
        sample = fields[1]

        sample_gene = sample +"_" + gene

        gene_counts[sample_gene] += 1


for item,count in gene_counts.items():
    if count == 1:

        homozygous.append(item)


    else:
        continue

for line in open(sv_calls):
        
    if line.startswith("#"):
        print line.strip()
    else:
        fields = line.strip().split("\t")
        
        gene = fields[0]
        sample = fields[1]

        sample_gene = sample +"_" + gene

        if sample_gene in homozygous:
            print line.strip()
            print line.strip()
        else:
            print line.strip()


