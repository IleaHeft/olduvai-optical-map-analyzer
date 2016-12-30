#! /usr/bin/env python

import pdb
import sys
from collections import defaultdict

directory=sys.argv[1]
sample = directory.split("/")[1]

shift_nick = int(sys.argv[2])

# either "HLS" or "CON1"
region_to_measure = sys.argv[3]

# File created in irys-analyzer.sh that gives reference nick positions of interest
duf_nicks=sys.argv[4]

# output folder
output_dir = sys.argv[5]

# Provide paths to each of these files
xmap=directory + "/alignmolvref_contig1.xmap"
r_cmap=directory + "/alignmolvref_contig1_r.cmap" 
q_cmap=directory + "/alignmolvref_contig1_q.cmap"




# Set up output file
if region_to_measure == "HLS":
    out_each_sample = open(sample + "-con2-con3-dist-plus" + str(shift_nick) + ".tab", mode = 'w')

    ref_dist_file = open(output_dir + "/ref-dist-hls.tab", mode = 'a')

else:
    out_each_sample = open(output_dir + "/" + sample + "-con1-region-dist-plus" + str(shift_nick) + ".tab", mode = 'w')

    ref_dist_file = open(output_dir + "/ref-dist-con1-region.tab", mode = 'a')

# Set up an output file for status updates/flags
status = open(output_dir + "/status.txt", mode = 'a')
#flagged_mol_file = open("nick-dist-flagged-mols.txt", mode = 'a')

#out_one_file = open("con2-con3-lengths-all-samples.tab", mode = 'w')

flagged_mols = defaultdict(tuple)
left_nicks = defaultdict(list)
right_nicks = defaultdict(list)
poor_align_right = []
poor_align_left = []
left_ref_qual = defaultdict(str)
right_ref_qual = defaultdict(str)
left_query_qual = defaultdict(str)
right_query_qual = defaultdict(str)
align_qual = defaultdict(str)
diff_bw_nick_nums = defaultdict(str)

#Step 1: Generate a dictionary of reference positions for nick sites of interest
gene_clades = defaultdict(list)
gene_nicks = defaultdict(list)
gene_strands = defaultdict(str)

for line in open(duf_nicks):
    fields = line.strip().split("\t")

    chrom = fields[0]
    nick_location = fields[1] + ".0"
    domain  = fields[6]
    clade = domain.split("_")[1]
    gene = fields[7]
    strand = fields[8]

    gene_clades[gene].append(gene)
    gene_strands[gene] = strand

    if "CON2" in domain or "CON3" in domain:
        if "CON2" in domain:
            gene_nicks[gene].append((clade,nick_location))
     
            gene_clades[gene].append(clade)
       
        else:
            if "CON3" not in gene_clades.get(gene):
                gene_clades[gene].append(clade)
    
                gene_nicks[gene].append((clade,nick_location))

    


## Print the reference distances to a files that can be used by the sv caller python script

#if region_to_measure == "HLS":
 #   for gene, nicks in gene_nicks.items():
  #      if len(nicks) == 2:
   #         nick1_ref = float(nicks[0][1])
    #        nick2_ref = float(nicks[1][1])
     #       ref_dist = abs(nick2_ref - nick1_ref)
      #      ref_dist_line = [sample, gene, str(nick1_ref), str(nick2_ref), str(ref_dist)]
            
       #     print >> ref_dist_file, "\t".join(ref_dist_line)
        

# Step 2: For every set of gene nicks, use the reference cmap file to identify the nick site number of the reference nicks of interst

for gene, nicks in gene_nicks.items():
    
    if len(nicks) != 2:
        if "CON2" in nicks[0]: 
            nick1_clade = nicks[0][0]
            nick1_ref = nicks[0][1]
            num_nicks_in_gene = len(nicks)
            print >> status, gene, "only one nick site"
        else:
            continue

    else:
        nick1_clade = nicks[0][0]
        nick1_ref = nicks[0][1]
     
        nick2_clade = nicks[1][0]
        nick2_ref = nicks[1][1]
        
        num_nicks_in_gene = len(nicks)

    for line in open(r_cmap):
        
        if line.startswith("#"):
            continue
        else:
            fields = line.strip().split("\t")
            query_id = fields[0]
            mol_length = fields[1]
            site_num = fields[3]
            ref_pos = fields[5]

            
            if ref_pos == nick1_ref:
                ref_nick1_num = int(site_num)
                
                
                if num_nicks_in_gene == 1:
                    if nick1_clade == "CON2":
                        ref_CON2_num = ref_nick1_num

                        if gene_strands.get(gene) == "+":
                            ref_CON2R1_num = ref_nick1_num + 1
                        else:
                            ref_CON2R1_num = ref_nick1_num - 1
                        
                        
                    if nick1_clade == "CON3":
                        ref_CON3_num = ref_nick1_num

                        if gene_strands.get(gene) == "+":
                            ref_CON3L1_num = ref_nick1_num - 1
                        else:
                            ref_CON3L1_num = ref_nick1_num + 1
                        
                
                # This is for calculating CON1 region distances - we take the reference nick number for CON2 and subtract one to get the nick number of one upstream - neg strand is accounted for by adding 1 to ref nick 2 num
                if nick1_clade == "CON2":
                    ref_CON2L1_num = ref_nick1_num - 1
                    ref_CON2_num = ref_nick1_num
                    
                    
                # This step shifts the nick site the specified number of nicks, if desired
                if nick1_clade == "CON3":
                    ref_nick1_num = int(ref_nick1_num - shift_nick)

            elif ref_pos == nick2_ref:
                ref_nick2_num = int(site_num)


                # This is for calculating CON1 region distances - if the nick2 clade is CON2, then the gene is on the negative strand, and to get the right nick for CON1 region measurement, we add 1 to the ref_nick2_num
                if nick2_clade == "CON2":
                    ref_CON2L1_num = ref_nick2_num +1
                    ref_CON2_num = ref_nick2_num

                # This step shifts the nick site the specified number of nicks, if desired
                if nick2_clade == "CON3":
                    if strand == "+":
                        ref_nick2_num = int(ref_nick2_num + shift_nick)
                    if strand == "-":
                        ref_nick2_num = int(ref_nick2_num - shift_nick)
            else:
                continue
    
    if region_to_measure == "HLS":
        if num_nicks_in_gene == 1:
            
            if nick1_clade == "CON2":
                if gene_strands.get(strand) == "+":
                    ref_nick1_num = ref_CON2_num
                    ref_nick2_num = ref_CON2R1_num
                else:
                    ref_nick1_num = ref_CON2R1_num
                    ref_nick2_num = ref_CON2_num

            
            if nick1_clade == "CON3":
                if gene_strands.get(strand) == "+":
                    ref_nick1_num = ref_CON3L1_num
                    ref_nick2_num = ref_CON3_num
                else:
                    ref_nick1_num = ref_CON3_num
                    ref_nick2_num = ref_CON3L1_num

    else:
        if nick1_clade == "CON2":
            ref_nick1_num = ref_CON2L1_num
            ref_nick2_num = ref_CON2_num
        else:
            ref_nick1_num = ref_CON2_num
            ref_nick2_num = ref_CON2L1_num

            
    # print reference distances
    for line in open(r_cmap):
        
        if line.startswith("#"):
            continue
        else:
            fields = line.strip().split("\t")
            query_id = fields[0]
            mol_length = fields[1]
            site_num = fields[3]
            ref_pos = fields[5]
            
            if int(site_num) == int(ref_nick1_num):
                ref_pos1 = ref_pos
            elif int(site_num) == int(ref_nick2_num):
                ref_pos2 = ref_pos
            else:
                continue

    ref_dist = abs(float(ref_pos2) - float(ref_pos1))
    print_ref_dist = [sample, gene, ref_pos1, ref_pos2, str(ref_dist)]
    print >> ref_dist_file, "\t".join(print_ref_dist)
    
    
    # Step 3: Look into xmap file alignment strings and identify molecules aligned to the reference nick of interest - as well as the relevant query nick sites in each molecule 


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
            strand = fields[7]
            align = fields[13]
            
            mol_conf[query_id] = confidence

            if confidence >= 8:
            # This strips the initial "(" and the ")" on the end, but looking at the string starting 1 position in from either end, then splits it up by the parentheses combination that separates each alignment pair
            # This is a list of comma separated pairs

                # Split the alignment string into a list of alignment pairs
                alignments = align[1:len(align)-1].split(")(")
                
                # Initialize empty lists
                ref_nicks_in_align_str = []
                query_nicks_in_align_str = []

                for align_pair in alignments:
                    align_pair = align_pair.split(",")
                    ref_nick_num = align_pair[0]
                    query_nick_num = align_pair[1]


                    ref_nicks_in_align_str.append(ref_nick_num)
                    query_nicks_in_align_str.append(query_nick_num)

                    if ref_nick_num == str(ref_nick1_num):
                        query_nick1_num = query_nick_num
                        nick_dict[query_id].append((query_nick1_num))
                        relevant_molecules.append(query_id)
                        #nick1_items.append((query_id, query_nick_num))

                    elif ref_nick_num == str(ref_nick2_num):
                        query_nick2_num = query_nick_num
                        
                        nick_dict[query_id].append((query_nick2_num))
                        
                        if query_id not in relevant_molecules:
                            relevant_molecules.append(query_id)
                        
                        #nick2_items.append((query_id, query_nick2_num))

                    else:
                        continue



                if query_id in relevant_molecules and len(nick_dict.get(query_id)) == 2:
 
                    # check that query nick 1 and query nick2 are adjacent to one another
                    
                    diff_bw_nicks =  abs(int(query_nick2_num) - int(query_nick1_num))
                    diff_bw_nick_nums[query_id] = str(diff_bw_nicks)
                    if diff_bw_nicks > 1:
                        status_line = ["FLAG, query nicks not adjacent", query_id, query_nick1_num, query_nick2_num] 
                        print >> status, "\t".join(status_line)
                        
                        flag_info = [sample, gene, query_id, query_nick1_num, query_nick2_num, str(diff_bw_nicks),str(confidence)]
                        flagged_mols[query_id] = "\t".join(flag_info)
                    

                    # Make dicionaries that hold all of the nicks to the left and the right of the reference nicks of interest
                    for ref_nick in ref_nicks_in_align_str:
                         if int(ref_nick) < int(ref_nick1_num):
                             
                             query_id_gene = query_id + "_" + gene

                             left_nicks[query_id_gene].append(ref_nick)
 
                         elif int(ref_nick) > int(ref_nick2_num):


                            query_id_gene = query_id + "_" + gene

                            right_nicks[query_id_gene].append(ref_nick)
                         
                         else:
                            continue

                    # Get the indices of the left and right reference nicks of interest (this could be done with list.index(), but doing it this way allows a check for duplicates of the reference nick in the alingment string)
                    
                    # Make empty lists to use in the check to make sure the nick doesn't occur more than once
                    index_ref1_check = []
                    index_ref2_check = []
                    
                    #for index, item in enumerate(ref_nicks_in_align_str):

                     #   if item == ref_nick1_num:
                      #      index_ref_nick1 = index
                            
                            # Part of built in check to make sure reference nick doesn't occur more than once
                       #     index_ref_nick1_check.append(index)

                       # elif item == ref_nick2_num:
                        #    index_ref_nick2 = index
                        
                            # Part of built in check, as above, but for reference nick two
                         #   index_ref_nick2_check.append(index)

                       # else:
                            #continue
                    
                    index_ref1 = ref_nicks_in_align_str.index(str(ref_nick1_num))
                    index_ref2 = ref_nicks_in_align_str.index(str(ref_nick2_num))

                    index_query1 = query_nicks_in_align_str.index(str(query_nick1_num))
                    index_query2 = query_nicks_in_align_str.index(str(query_nick2_num))
                    
                    # Check that reference nicks aligned to left and right of nicks of interest are one value larger and smaller
                    
                    

                      # Check that reference nick to the left of nick1 is one value smaller 
                      # Only run this check for molecules that have at least one nick to the left
                    if query_id_gene in left_nicks.keys():
                        
                        if strand == "-":
                            adj = -1
                        else:
                            adj = 1

                        # check the order of the query nicks to the left
                        if int(query_nicks_in_align_str[index_query1 - 1]) == int(query_nick1_num) - adj:
                            
                            left_query_qual = "G"

                            status_line = [sample,query_id, strand, str(query_nick1_num), str(query_nick2_num), "query left nick minus 1: PASS", align]
                            print >> status, "\t".join(status_line)
                            
                        else:
                            left_query_qual = "P"

                            status_line = [sample,query_id, strand, str(query_nick1_num), str(query_nick2_num), "query left nick minus 1: FAIL", align]
                            print >> status, "\t".join(status_line)


                            poor_align_left.append(query_id)




                        # check the order of the reference nicks to the left
                        if int(ref_nicks_in_align_str[index_ref1 - 1]) == int(ref_nick1_num) - 1:
                            
                            left_ref_qual= "G"

                            status_line = [sample,query_id, strand,str(ref_nick1_num),str(ref_nick2_num),"ref left nick minus 1: PASS", align]
                            print >> status, "\t".join(status_line)
                        
                        else:
                            status_line = [sample,query_id, strand,str(ref_nick1_num),str(ref_nick2_num), "ref left nick minus 1: FAIL", align]
                            print >> status, "\t".join(status_line)
                            
                            left_ref_qual = "P"
                        # Check the possibility of false negative on molecule giving appearance of skip on reference
                        if left_query_qual == "G" and left_ref_qual == "P" and len(left_nicks.get(query_id_gene)) > 1: #only bother checking if the query alignment nicks to the left is in order
                            # the frist condition checks the 2 positions to the left of ref_nick1_num, second condition checks the position in the query 2 positions to the left of query_nick1_num

                            if strand == "-":
                                adj2 = -2
                            else:
                                adj2 = 2

                            if int(ref_nicks_in_align_str[index_ref1 - 2]) == int(ref_nick1_num) - 3 and int(query_nicks_in_align_str[index_query1 - 2]) == int(query_nick1_num) - adj2:
                                
                                left_ref_qual = "S"

                            else:
                                left_ref_qual = "P"


                        # Check the possibility of a false positive or novel nick in query

                        if left_ref_qual == "G" and left_query_qual == "P" and len(left_nicks.get(query_id_gene)) > 1: #only bother checking if the query alignment nicks to the left is in order
                            # the frist condition checks the 2 positions to the left of ref_nick1_num, second condition checks the position in the query 2 positions to the left of query_nick1_num

                            if strand == "-":
                                adj3 = -3
                            else:
                                adj3 = 3

                            if int(query_nicks_in_align_str[index_query1 - 2]) == int(query_nick1_num) - adj3 and int(ref_nicks_in_align_str[index_ref1 - 2]) == int(ref_nick1_num) - 2:
                                
                                left_query_qual = "S"

                            else:
                                left_query_qual = "P"


                        poor_align_left.append(query_id)
                    else:
                        left_ref_qual = "U"
                        left_query_qual = "U"


                      # Check that reference nick to the right of nick2 is one value larger
                      # Only run this check for molecules that have at least one nick to the right
                    if query_id_gene in right_nicks.keys():

                        if strand == "-":
                            adj = -1
                        else:
                            adj = 1

                        # Check the query alignment
                        
                        
                        
                        if int(query_nicks_in_align_str[index_query2 + 1]) == int(query_nick2_num) + adj:

                            right_query_qual = "G"

                            status_line = [query_id, strand, str(query_nick1_num), str(query_nick2_num), "query right nick plus 1: PASS", align]
                            print >> status, "\t".join(status_line)
                        else:

                            right_query_qual = "P"
                            status_line = [query_id, strand, str(query_nick1_num), str(query_nick2_num), "query right nick plus 1: FAIL", align]
                            print >> status, "\t".join(status_line)


                            poor_align_right.append(query_id)




                        # check the reference alignmnent
                        if int(ref_nicks_in_align_str[index_ref2 + 1]) == int(ref_nick2_num) + 1:
                            
                            right_ref_qual = "G"

                            status_line = [query_id, strand,str(ref_nick1_num),str(ref_nick2_num), "ref right nick plus 1: PASS", align]
                            print >> status, "\t".join(status_line)
                            
                        else:
                            right_ref_qual = "P"

                            status_line = [query_id, strand,str(ref_nick1_num),str(ref_nick2_num), "ref right nick plus 1: FAIL", align]
                            print >> status, "\t".join(status_line)

                            poor_align_right.append(query_id)

                        # Check the possibility of false negative on molecule giving appearance of skip on reference
                        if right_query_qual == "G" and right_ref_qual == "P" and len(right_nicks.get(query_id_gene)) > 1: #only bother checking if the query alignment nicks to the left is in order
                            # the frist condition checks the 2 positions to the left of ref_nick1_num, second condition checks the position in the query 2 positions to the left of query_nick1_num

                            if strand == "-":
                                adj2 = -2
                            else:
                                adj2 = 2

                            if int(ref_nicks_in_align_str[index_ref2 + 2]) == int(ref_nick2_num) + 3 and int(query_nicks_in_align_str[index_query2 + 2]) == int(query_nick2_num) + adj2:
                                
                                right_ref_qual = "S"

                            else:
                                right_ref_qual = "P"

                        # Check the possibility of false positive or novel nick on query 
                        if right_query_qual == "P" and right_ref_qual == "G" and len(right_nicks.get(query_id_gene)) > 1: #only bother checking if the query alignment nicks to the left is in order
                            # the frist condition checks the 2 positions tio the left of ref_nick1_num, second condition checks the position in the query 2 positions to the left of query_nick1_num

                            if strand == "-":
                                adj3 = -3
                            else:
                                adj3 = 3

                            if int(query_nicks_in_align_str[index_query2 + 2]) == int(query_nick2_num) + adj3 and int(ref_nicks_in_align_str[index_ref2 + 2]) == int(ref_nick2_num) + 2:
                                
                                right_query_qual = "S"

                            else:
                                right_query_qual = "P"




                        poor_align_left.append(query_id)

                    else:
                        right_ref_qual = "U"
                        right_query_qual = "U"

                    align_qual[query_id] = left_ref_qual + left_query_qual + "-" + right_ref_qual + right_query_qual 
                     # Do the check to make sure neither of the reference nicks occurs more than once
                    if len(index_ref1_check) > 1:
                        print >> status, query_id, "FLAG: Reference nick1 present more than once in alignment string"
                    elif len(index_ref2_check) > 1:
                        print >> status, query_id, "FLAG: Reference nick2 present more than once in alignment string"
                    else:
                        continue
                        



    query_pos_dict = defaultdict(list)
    
    for line in open(q_cmap):
        if line.startswith("#"):
            continue
        else:
            fields = line.strip().split("\t")

            query_id = fields[0]
            mol_length = fields[1]
            site_num = fields[3]
            position = fields[5]

            if query_id in nick_dict.keys():
                if len(nick_dict.get(query_id)) == 2:
                    nick_nums = nick_dict.get(query_id)
                    positions = [] 
                    if site_num in nick_nums:

                        if site_num == nick_nums[0]:
                            position_nick1 = position
                            query_pos_dict[query_id].append(position_nick1)
                            positions.append(position_nick1)
                        elif site_num == nick_nums[1]:
                            position_nick2 = position
                            query_pos_dict[query_id].append(position_nick2)
                            positions.append(position_nick2)
                        else:
                            continue
    
    for query_id, positions in query_pos_dict.items():
        

        dist = abs(float(positions[0]) - float(positions[1]))
        confidence = str(mol_conf.get(query_id))
      
        query_id_gene = query_id + "_" + gene   

        if query_id_gene in left_nicks.keys():

            num_left_nicks = len(left_nicks.get(query_id_gene))
        else:
            num_left_nicks = 0
            
            # set align_qual to "U" for NA because, since there aren't nicks to the left, we can't check alignment of nicks to the left (same for right nick conditions)
            align_qual_left = "U"

        if query_id_gene in right_nicks.keys():
            num_right_nicks = len(right_nicks.get(query_id_gene))
        else:
            num_right_nicks = 0
            
            align_qual_right = "U"

        # set the align_qual string to be printed to be the combination of the left alignment and right alignment status

        #align_qual = align_qual_left + align_qual_right
        if len(gene_nicks.get(gene)) == 2:
            nicks_used = "std"
        else:
            nicks_used = "alt"


        align_qual_flags = align_qual.get(query_id)

        to_print_one_file = [sample,gene,query_id,confidence,diff_bw_nick_nums.get(query_id),align_qual_flags,str(num_left_nicks),str(num_right_nicks),str(dist),nicks_used]
        print "\t".join(to_print_one_file)
    
        #print >> out_each_sample, "\t".join(to_print_one_file)

            
        #print >> flagged_mol_file, flagged_mols.get(query_id)
