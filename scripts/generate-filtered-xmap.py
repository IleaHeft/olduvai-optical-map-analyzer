#! /usr/bin/env python

import sys
import pdb
from collections import defaultdict

sample_dir = sys.argv[1]
sample = sys.argv[2]
xmap = sys.argv[3]
conf_spread = int(sys.argv[4])


mol_conf = defaultdict(list)
to_filter_out = []
to_keep = defaultdict(float)

# Set up output files
filtered_out = open(sample_dir + "/" + sample + "-mols-with-ambiguous-alignment.txt",mode = 'w')
filtered_xmap = open(sample_dir + "/" + sample + "-filt.xmap",mode = 'w')


for line in open(xmap):
    if line.startswith("#"):
        continue
    else:
        fields = line.strip().split("\t")
        query_id = fields[1]
        confidence = float(fields[8])
        strand = fields[7]
        align = fields[13]

        mol_conf[query_id].append(confidence)

for query_id, conf_scores in mol_conf.items():
    conf_scores = sorted(conf_scores)
 
    if len(conf_scores) > 1:
        max_conf = max(conf_scores)

        for score in conf_scores:
            if score != max_conf:
                diff_from_max = float(max_conf) - float(score)

                if diff_from_max <= conf_spread:
                    to_filter_out.append(query_id)
                    print >> filtered_out, sample,"\t",query_id,"\t",conf_scores
                else:
                    to_keep[query_id] = max_conf

        # Filter out molecules where the identical maximum score occurs more than once
        num_max_scores = 0
        for score in conf_scores:
            if score == max_conf:
                num_max_scores += 1
        if num_max_scores > 1:
            to_filter_out.append(query_id)
            print >> filtered_out, sample,"\t",query_id,"\t",conf_scores
                
    else:
        to_keep[query_id] = conf_scores[0]

for line in open(xmap):
    if line.startswith("#"):
        print >> filtered_xmap, line.strip()
    else:
        fields = line.strip().split("\t")
        query_id = fields[1]
        confidence = float(fields[8])
        strand = fields[7]
        align = fields[13]

        if query_id in to_filter_out:
            continue
        elif confidence == to_keep[query_id]:
            print >> filtered_xmap, line.strip()
        else:
            continue


