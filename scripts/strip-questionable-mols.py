#! /usr/bin/env python

import sys
import pdb

mols_to_strip = sys.argv[1]
results_to_clean = sys.argv[2]

mols_to_strip_list = []

for line in open(mols_to_strip):
    fields = line.strip().split("\t")
    sample = fields[0].strip()
    mol = fields[1].strip()

    sample_mol = sample + "_" + mol
    mols_to_strip_list.append(sample_mol)


for line in open(results_to_clean):
    
    fields = line.strip().split("\t")
    sample = fields[0]
    mol = fields[2]

    sample_mol = sample + "_" + mol
    
    if sample_mol not in mols_to_strip_list:
        print "\t".join(fields)

    else:
        continue

