#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import pandas as pd
import glob
import re

parser = argparse.ArgumentParser(description="Derive the portion of each chimeric gene actually aligned in each orthtogroup it belongs to")
parser.add_argument("--input_dir", "-i", required=True, metavar="input_dir", help="Directory containing all the multiple alignments between the chimeric genes and the respective outgroups")
parser.add_argument("--alignment", "-a", required=True, metavar="alignment", help="List of alignment files prefix in the format chimeric_geneID-orthogroupID")
parser.add_argument("--exon_positions_dir", "-e", required=True, metavar="exon_positions_dir", help="Directory containing all species files with the aminoacidic positions corresponding to each exon")
#parser.add_argument("--chimeric_gene", "-g", required=True, metavar="chimeric_gene", help="GeneID of the chimeric gene in the alignment")
#parser.add_argument("--orthogroup_id", "-og", required=True, metavar="orthogroup_id", help="OrthogroupID of the orthogroup to which the chimeric gene belongs")
parser.add_argument("--overlap_stringency", "-s", required=True, metavar="overlap_stringency", help="Minumum percentage of aligned genes in OG required to define the overlap")
parser.add_argument("--output", "-o", required=True, metavar="output_dir", help="path to output file")

#Read arguments
args = parser.parse_args()
my_input_dir = args.input_dir
my_exon_positions_dir = args.exon_positions_dir
my_overlap_stringency = float(args.overlap_stringency)
my_output_file = args.output
with open(args.alignment, "r") as alignment_prefix_file:
  my_alignment_prefix_list = alignment_prefix_file.readline().rstrip().split(";") #list of chimeric_geneID;orthogroupID

#Define functions
#This function gets the number of sligned position in the string to which the tested AA belongs to.
def aligned_pos_in_string(aln_list, x):
  num = 0
  pos_list_left = list(range(0,x))[::-1] #the right boundary is excluded. rever the list
  pos_list_right = list(range(x+1,len(aln_list)-1))
  for pos in pos_list_left:
    if aln_list[pos] != "-":
      num = num+1
    else:
      break
  for pos in pos_list_right:
    if aln_list[pos] != "-":
      num = num+1
    else:
      break
  return(num+1) #the +1 counts the position itself

#Generate dataframe with all the correspondence between exon num and aa_coords for all genes and all species.
all_exon_pos_files = glob.glob(my_exon_positions_dir+"/*_refprots_exons_pos.tab")
all_exon_pos_df = pd.DataFrame()
for my_file in all_exon_pos_files:
  my_df = pd.read_table(my_file, sep="\t", index_col=False, header=None, names=["geneID", "exon_num","aa_coords", "genomic_coords", "strand"]) #NB: the GeneID col is still in the format ProtID|GeneID
  all_exon_pos_df = pd.concat([all_exon_pos_df, my_df])
all_exon_pos_df["geneID"] = [re.sub(".*\|", "", str(element)) for element in list(all_exon_pos_df["geneID"])] #modify geneID column

#print header in output fiels
with open(my_output_file, "a") as my_output:
  my_output.write("chimeric_geneID\torthogroup_id\tfirst_aligned_aa\tlast_aligned_aa\tfirst_aligned_ex\tlast_aligned_ex\n")

#cicle on each multiple alignment
for my_alignment_prefix in my_alignment_prefix_list:
  my_chimeric_gene = my_alignment_prefix.split("-")[0]
  my_orthogroupID = my_alignment_prefix.split("-")[1]
  my_alignment_file = my_input_dir+"/"+my_alignment_prefix+"-multiple_aln"
  #print chimeric gene (just to debug)
  print(my_chimeric_gene+"\t"+my_orthogroupID)

  #Get a list of all the entries in the alignment 
  my_aln = list(list(SeqIO.parse(my_alignment_file, "fasta")))
  total_genes = len(my_aln) #Count total amount of genes in the orthogroup
  chimeric_gene_aln_seq = [element for element in my_aln if element.id == my_chimeric_gene][0] #Isolate the chimeric gene entry

  #Create a dataframe with row = alignment_position, column = gene
  alignment_df = pd.DataFrame()
  for entry in my_aln:
    alignment_df[entry.id] = list(entry.seq)
  transposed_alignment_df = alignment_df.transpose() #Transpose the dataframe (row=gene, column=alignment position)

  #Get a dictionary with key=alignment position, value=percentage of genes that align there (i.e. position != "-")
  max_aln_positions = transposed_alignment_df.shape[0] #number of species
  percent_aln_positions = transposed_alignment_df.apply(lambda x: len([element for element in x if element != "-"])/max_aln_positions, axis=0) #this returns a series
  position_percent_aln_dict = percent_aln_positions.to_dict()
  #Get a dictionary with key=alignment position, value=length of the aligned  fragment to which the position belongs
  aligned_pos_len_list = [aligned_pos_in_string("".join(list(entry.seq)), x) for x in list(range(len("".join(list(entry.seq)))))]
  aligned_pos_num_dict = {index: value for index,value in enumerate(aligned_pos_len_list)}
  
  #Select the lowest and highest key which are aligned and in chimeric gene for which value >= %stringency
  #first and last aligned aminoacid, respectively
  chimeric_gene_aln_pos = transposed_alignment_df.loc[my_chimeric_gene,]
  chimeric_gene_aligned_pos = list(chimeric_gene_aln_pos[chimeric_gene_aln_pos != "-"].index.values)
  valid_positions = [key for key in list(position_percent_aln_dict.keys()) if key in chimeric_gene_aligned_pos]
  valid_positions_with_length = [position for position in valid_positions if aligned_pos_num_dict[position] >= 10] #the position is in a string of at least n aligned AA
  valid_positions_with_coverage = [position for position in valid_positions_with_length if position_percent_aln_dict[position] >= my_overlap_stringency] #at least one positions that aligns with a certain percentage of the orthogroup
  if len(valid_positions_with_coverage) >=1:
    first_aligned_pos = min(valid_positions_with_coverage)
    last_aligned_pos = max(valid_positions_with_coverage)

    #Translate from alignment position to chimeric gene position. Create a dictionary between the two "coordinate" systems (the chimeric gene aln with and without indels)
    aln_chimeric_gene_pos_dict = {key: value for value, key in enumerate(valid_positions)}
    first_aligned_aa = aln_chimeric_gene_pos_dict[first_aligned_pos]+1 #first aa in the fragment. Add one because the other aa_pos are in a 1-based system
    last_aligned_aa = aln_chimeric_gene_pos_dict[last_aligned_pos]+1 #last aa in the fragment. Add one because other aa_post are in a 1-based system

    #Subset the exon positions dataframe to only the chimeric gene of interest
    chimeric_gene_exon_pos = all_exon_pos_df[all_exon_pos_df["geneID"]==my_chimeric_gene]
    chimeric_gene_exon_pos_dict = pd.Series(chimeric_gene_exon_pos.exon_num.values, index=chimeric_gene_exon_pos.aa_coords).to_dict() #dictionary with key=aa_pos, value=exon_num
    chimeric_gene_exon_pos_ranges = [range(int(element.split("-")[0]), int(element.split("-")[1])) for element in list(chimeric_gene_exon_pos_dict.keys())]
    overlapping_exons = [element for element in chimeric_gene_exon_pos_ranges if set(range(first_aligned_aa, last_aligned_aa)).intersection(element)]
    if len(overlapping_exons) >= 1:
      first_overlapping_exon = overlapping_exons[0]
      last_overlapping_exon = overlapping_exons[len(overlapping_exons)-1]
      #Transform ranges back to coordinates
      first_overlapping_exon_aa = str(first_overlapping_exon[0])+"-"+str(first_overlapping_exon[len(first_overlapping_exon)-1]+1)
      last_overlapping_exon_aa = str(last_overlapping_exon[0])+"-"+str(last_overlapping_exon[len(last_overlapping_exon)-1]+1)
      #Intersect with the aminoacid positions to get the correspondent exon numer
      first_overlapping_exon_num = chimeric_gene_exon_pos_dict[first_overlapping_exon_aa]
      last_overlapping_exon_num = chimeric_gene_exon_pos_dict[last_overlapping_exon_aa]
    else:
      first_overlapping_exon_aa = "NA"; last_overlapping_exon_aa = "NA"; first_overlapping_exon_num = "NA"; last_overlapping_exon_num = "NA"
  else:
    first_aligned_aa = "NA"; last_aligned_aa = "NA"; first_overlapping_exon_aa = "NA"; last_overlapping_exon_aa = "NA"; first_overlapping_exon_num = "NA"; last_overlapping_exon_num = "NA"

  #Save output file with header:
  #1. Chimeric geneID
  #2. OG ID
  #3. First aligned aa
  #4. Last aligned aa
  #5. First aligned exon (exon num: aa coords)
  #6. Last aligned exon (exon num: aa coords)
  final_string = "%s\t%s\t%s\t%s\t%s:%s\t%s:%s\n" % (my_chimeric_gene, my_orthogroupID, first_aligned_aa, last_aligned_aa, first_overlapping_exon_num, first_overlapping_exon_aa, last_overlapping_exon_num, last_overlapping_exon_aa)
  with open(my_output_file, "a") as my_output:
     my_output.write(final_string)

  print(final_string)
