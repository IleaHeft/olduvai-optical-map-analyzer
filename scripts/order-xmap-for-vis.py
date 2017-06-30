#! /usr/bin/env bash

import sys
import pdb
from collections import defaultdict

xmap = sys.argv[1]
order = sys.argv[2]

mol_dict = defaultdict(int)

for index,line in enumerate(open(order)):
    fields = line.strip().split("\t")
    mol_id = fields[2]

    num = index + 1
    
    mol_dict[mol_id] = num
    

for line in open(xmap):
    if line.startswith("#"):
        print line.strip()

    else:
        fields = line.strip().split("\t")
        ref_start = fields[5]
        mol_id = fields[1]
        
        if mol_id in mol_dict.keys():

            new = mol_dict.get(mol_id)
        
            fields[5] = str(new)
            
            to_print = "\t".join(fields)
            print to_print
