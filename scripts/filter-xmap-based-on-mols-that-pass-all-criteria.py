#! /usr/bin/env python

import sys
import pdb
from collections import defaultdict

mols_passing = sys.argv[1]
xmap = sys.argv[2]

mols_to_retain = []
mols_printed = []

# Set up output files

for line in open(mols_passing):
    mols_to_retain.append(line.strip())



for line in open(xmap):
    if line.startswith("#"):
        print line.strip()
    else:
        fields = line.strip().split("\t")
        query_id = fields[1]
        
        if query_id in mols_to_retain:
           print line.strip()

           mols_printed.append(query_id)
        else:
            continue

