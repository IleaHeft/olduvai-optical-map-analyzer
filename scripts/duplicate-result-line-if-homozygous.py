#! /usr/bin/env python

import sys
import pdb

zygosity_calls = sys.argv[1]

for line in open(zygosity_calls):
    if line.startswith("#"):
        print line.strip()
    else:
        fields = line.strip().split("\t")
        num_calls = fields[3]
        if num_calls == "1":
            print line.strip()
            print line.strip()
        else:
            print line.strip()
