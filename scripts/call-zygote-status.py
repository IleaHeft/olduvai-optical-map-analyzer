#! /usr/bin/env python

import sys
import pdb

sv_calls = sys.argv[1]


for line in open(sv_calls):
    fields = line.strip().split("\t")
    sample = fields[0]
    gene = fields[1]
    count = fields[2]
    types = fields[3]
    sizes = fields[4]

    types_list = types.split(",")
    
    if "nd" in types:

        to_print = [sample, gene, "nd", "na", "na", "na"]
        print "\t".join(to_print)        
    
    else:
        if "ins" not in types and "del" not in types and "nd" not in types:
            genotype = "homozygous-reference"

        elif count == "1" and "ref" not in types:
            genotype = "homozygous-" + types
        elif count == "2":
            if "ref" in types:
      
                for item in types_list:
                    if item == "ref":
                        type1 = "ref"
                    else:
                        type2 = item
            elif "ins" in types and "del" not in types:
                type1 = "ins"
                type2 = "ins"

            elif "del" in types and "ins" not in types:
                type1 = "del"
                type2 = "del"

            else:
                type1 = "ins"
                type2 = "del"

            genotype = "heterozygous-"+ type1 + "," + type2 
        elif count == "3":
            
            ref_count = 0
            
            for item in types_list:
                if item == "ref":
                    ref_count += 1
                else:
                    other_type = item
            if ref_count == 2:

                genotype = "heterozygous-ref," + other_type 
            else:
                genotype = "heterozygous-" + types
        
        else: 
            genotype = "more than two"



        to_print = [sample, gene, genotype, count, types, sizes]
        print "\t".join(to_print)        

