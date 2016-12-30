#! /usr/bin/env python

import pdb
import sys
from collections import defaultdict

directory=sys.argv[1]
sample = directory.split("/")[1]

#shift_nick = int(sys.argv[2])
# Provide paths to each of these files
xmap=directory + "/alignmolvref_contig1.xmap"
r_cmap= directory + "/alignmolvref_contig1_r.cmap" 
q_cmap=directory + "/alignmolvref_contig1_q.cmap"


# Give known reference position of nicks of interest
con2_con3_nicks="con2-con3-nick-sites.bed"


# Set up output file
#out_each_sample = open(sample + "relevant-molecules.tab", mode = 'w')

#out_one_file = open("con2-con3-lengths-all-samples.tab", mode = 'w')

#Generate a dictionary of reference positions for nick sites of interest
gene_nicks = defaultdict(list)
for line in open(con2_con3_nicks):
    fields = line.strip().split("\t")

    chrom = fields[0]
    nick_location = fields[1] + ".0"
    domain  = fields[6]
    clade = domain.split("_")[1]
    gene = fields[7]
    strand = fields[8]

    if "CON2" in domain or "CON3" in domain:
        gene_nicks[gene].append((clade,nick_location))


# Print cmap header

for line in open(q_cmap):
    if line.startswith("#"):
        print line.strip()
    else:
        continue




# For every set of gene nicks, use the reference cmap file to identify the nick site number of the reference nicks of interst
for gene, nicks in gene_nicks.items():
    nick1_clade = nicks[0][0]
    nick1_ref = nicks[0][1]
    
    nick2_clade = nicks[1][0]
    nick2_ref = nicks[1][1]

#    for line in open(r_cmap):
        
 #       if line.startswith("#"):
  #          continue
   #     else:
    #        fields = line.strip().split("\t")
     #       query_id = fields[0]
      #      mol_length = fields[1]
       #     site_num = fields[3]
        #    ref_pos = fields[5]

         #   if ref_pos == nick1_ref:
          #      ref_nick1_num = int(site_num)
           #     if nick1_clade == "CON3":
            #        ref_nick1_num = ref_nick1_num - shift_nick
                #    if strand == "+":
                        
                 #       ref_nick1_num = ref_nick1_num + shift_nick
                    
                  #  if strand == "-":
                   #     ref_nick1_num = ref_nick1_num - shift_nick

          #  elif ref_pos == nick2_ref:
           #     ref_nick2_num = int(site_num)

            #    if nick2_clade == "CON3":
             #       if strand == "+":
              #          ref_nick2_num = ref_nick2_num + shift_nick
               #     if strand == "-":
                #        ref_nick2_num = ref_nick2_num - shift_nick
           # else:
            #    continue
    
    # Initialize dictionary: 
    ## keys = molecule IDs for molecules with nick aligned to reference nick of interest,
    ## values, the nick site number in the molecule that aligns with the reference nick of interest
    
    nick_dict=defaultdict(list)
    
    # Initialize empty list that will hold the molecule IDs of molecules aligned to the reference nick of interest
    relevant_molecules = []
    
    # Initialize dictionary:
    ## keys = molecule IDs for molecules with nick aligned to reference nick of interest,
    ## values = confidence values for the alignment of that molecule

    mol_conf=defaultdict(float)

    for line in open(xmap):
        if line.startswith("#"):
            continue
        else:
            fields = line.strip().split("\t")
            
            query_id = fields[1]
            confidence = float(fields[8])
            ref_align_start = float(fields[5])
            ref_align_end = float(fields[6])
            
            align = fields[13]
            
            mol_conf[query_id] = confidence

            if confidence >= 8:
            # This strips the initial "(" and the ")" on the end, but looking at the string starting 1 position in from either end, then splits it up by the parentheses combination that separates each alignment pair
            # This is a list of comma separated pairs
                
                # molecule spans both nick 1 and nick 2
                if ref_align_start <= float(nick1_ref) and ref_align_end >= float(nick2_ref):
                    relevant_molecules.append(query_id)

                # Molecule spans nick 1, but doesn't have to span nick 2
                elif ref_align_start <= float(nick1_ref) and ref_align_end >= float(nick1_ref):
                 
                    relevant_molecules.append(query_id)
                # Molecule spans nick 2, but doesn't have to span nick 1

                elif ref_align_start <= float(nick2_ref) and ref_align_end >= float(nick2_ref):

                    relevant_molecules.append(query_id)
                
                else:
                    continue
                 
                 #print query_id, ref_align_start, ref_align_end,gene,nick1_ref, nick2_ref
#                alignments = align[1:len(align)-1].split(")(")
#                
#                for align_pair in alignments:
#                    align_pair = align_pair.split(",")
#                    ref_nick_num = align_pair[0]
#                    query_nick_num = align_pair[1]
#
#                    if ref_nick_num == str(ref_nick1_num):
#                        query_nick1_num = query_nick_num
#                        nick_dict[query_id].append((query_nick1_num))
#                        relevant_molecules.append(query_id)
#                        #nick1_items.append((query_id, query_nick_num))
#                    
#                    elif ref_nick_num == str(ref_nick2_num):
#                        query_nick2_num = query_nick_num
#                        
#                        nick_dict[query_id].append((query_nick2_num))
#                        
#                        if query_id not in relevant_molecules:
#                            relevant_molecules.append(query_id)
#                        
#                        #nick2_items.append((query_id, query_nick2_num))
#
#                    else:
#                        continue
#    
#    query_pos_dict = defaultdict(list)
#    
    for line in open(q_cmap):
        if line.startswith("#"):
            continue
        else:
            fields = line.strip().split("\t")

            query_id = fields[0]
            mol_length = fields[1]
            site_num = fields[3]
            position = fields[5]

            if query_id in relevant_molecules:
                print line.strip()
           
           
#           if query_id in nick_dict.keys():
#                if len(nick_dict.get(query_id)) == 2:
#                    nick_nums = nick_dict.get(query_id)
#                    positions = [] 
#                    if site_num in nick_nums:
#
#                        if site_num == nick_nums[0]:
#                            position_nick1 = position
#                            query_pos_dict[query_id].append(position_nick1)
#                            positions.append(position_nick1)
#                        elif site_num == nick_nums[1]:
#                            position_nick2 = position
#                            query_pos_dict[query_id].append(position_nick2)
#                            positions.append(position_nick2)
#                        else:
#                            continue
#    
#    for query_id, positions in query_pos_dict.items():
#        length = abs(float(positions[0]) - float(positions[1]))
#        confidence = str(mol_conf.get(query_id))
#        to_print= [gene,query_id,confidence,str(length)]
#        print >> out_each_sample, "\t".join(to_print)
#        
#        to_print_one_file = [sample,gene,query_id,confidence,str(length)]
#        print "\t".join(to_print_one_file)
#    # Need to filter by confidence (start with a confidence level of 8)        
