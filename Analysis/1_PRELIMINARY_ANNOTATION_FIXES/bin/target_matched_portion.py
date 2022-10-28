#!/usr/bin/env python3

import sys
#NB: this has to be run in my Broccoli env to work. Otherwise, decomment the following line (which is ugly)
#sys.path.extend(['/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python36.zip', '/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python3.6', '/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python3.6/lib-dynload', '/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python3.6/site-packages'])

import argparse
from Bio import SeqIO
import re

parser = argparse.ArgumentParser(description="Compute the percentage of overlap (i.e. not gapped positions) between two aligned proteins")
parser.add_argument("--alignment", "-i", required=True, metavar="alignment_file", help="fasta file with the alignes proteins")

#NB: it is important that the second protein in the alignment is the common one between the two inputs.

def matched_portion(my_aln):
  query_aln_seq = my_aln[0].seq
  target_aln_seq = my_aln[1].seq
  query_matched_pos = [pos for pos,char in enumerate(query_aln_seq) if char!="-" and char!="X"]
  #this X screws things up: is it an ambigous aminoacid?

  matched_intervals = []
  original_start = query_matched_pos[0]
  start = query_matched_pos[0]
  count = 1
  for element in query_matched_pos:
    if element-start > 15:
      stop = query_matched_pos[query_matched_pos.index(element)-1] #get the preceding element in the list
      matched_intervals.append((original_start, stop))
      original_start = element
    else:
      #if element == query_matched_pos[-1]:
      if count == len(query_matched_pos):
        matched_intervals.append((original_start, element))
    start=element
    count=count+1

  matched_intervals_lengths = [element[1]-element[0] for element in matched_intervals] #length of the gaps
  #ordered_intervals_index = sorted(range(len(matched_intervals_lengths)), key=lambda k: matched_intervals_lengths[k], reverse=True)
  #ordered_matched_intervals = [x for _,x in sorted(zip(ordered_intervals_index,matched_intervals))]
  #chosen_interval = ordered_matched_intervals[0]
  
  #I need an alternative method because that sorting was screwing things up.
  my_max = max(matched_intervals_lengths) 
  length_dict = {element:element[1]-element[0] for element in matched_intervals}
  chosen_interval = [key for key in list(length_dict.keys()) if length_dict[key]==my_max][0]

  #select the sequence corresponding to these indexes in the target alignment and remove the gaps
  human_matched_seq = "".join([element for element in target_aln_seq[chosen_interval[0]:chosen_interval[1]+1] if element!="-"])
  #select the position in the original sequence which are matched 
  original_human_seq = "".join([element for element in target_aln_seq if element!="-"])
  modified_human_seq = re.sub(human_matched_seq, ["_"*len(human_matched_seq)][0], original_human_seq)
  #get the indexes of the positions matching "_"  
  matched_human_pos = [pos for pos,char in enumerate(modified_human_seq) if char=="_"]

  if len(matched_human_pos) != 0:
    matched_human_interval = (min(matched_human_pos)+1, max(matched_human_pos)+1) #tuple
  else:
    matched_human_interval=(0,0) #this is because Spu is stupid

  interval_to_write = str(matched_human_interval[0])+":"+str(matched_human_interval[1]) #transform into stringbeginning:end
  human_protein_length = len(original_human_seq) #derive protein length
  #write all information to output
  print("%s\t%s\t%d\t%d\t%s\t%d" % (my_aln[0].id, my_aln[1].id, matched_human_interval[0], matched_human_interval[1], interval_to_write, human_protein_length))
 #this prints to standard output

def main():
  args = parser.parse_args() #read arguments
  my_aln = list(list(SeqIO.parse(args.alignment, "fasta")))
  matched_portion(my_aln)

main()
