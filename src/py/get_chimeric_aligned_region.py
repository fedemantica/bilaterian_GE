#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import pandas as pd
import glob
import re

parser = argparse.ArgumentParser(description="Derive the portion of each chimeric gene actually aligned in each orthtogroup it belongs to")
parser.add_argument("--alignment", "-a", required=True, metavar="alignment", help="Multiple alignment among all the genes in the orthogroup containing chimeric genes")
parser.add_argument("--exon_positions_dir", "-e", required=True, metavar="exon_positions_dir", help="File with the aminoacidic positions corresponding to each exon")
parser.add_argument("--chimeric_gene", "-g", required=True, metavar="chimeric_gene", help="GeneID of the chimeric gene in the alignment")
parser.add_argument("--orthogroup_id", "-og", required=True, metavar="orthogroup_id", help="OrthogroupID of the orthogroup to which the chimeric gene belongs")
parser.add_argument("--overlap_stringency", "-s", required=True, metavar="overlap_stringency", help="Minumum percentage of aligned genes in OG required to define the overlap")
#parser.add_argument("--output", "-o", required=True, metavar="output_dir", help="path to output file")

#Read arguments
args = parser.parse_args()
my_alignment_file = args.alignment
my_exon_positions_dir = args.exon_positions_dir
my_chimeric_gene = args.chimeric_gene
my_orthogroupID = args.orthogroup_id
my_overlap_stringency = float(args.overlap_stringency)
#my_output = args.output

#Get a list of all the entries in the alignment 
my_aln = list(list(SeqIO.parse(my_alignment_file, "fasta")))
total_genes = len(my_aln) #Count total amount of genes in the orthogroup
chimeric_gene_aln_seq = [element for element in my_aln if element.id == my_chimeric_gene][0] #Isolate the chimeric gene entry

#Create a dataframe with row = alignment_position, column = gene
alignment_df = pd.DataFrame()
for entry in my_aln:
  alignment_df[entry.id] = list(entry.seq)
transposed_alignment_df = alignment_df.transpose() #Transpose the dataframe

#Get a dictionary with key=alignment position, value=percentage of genes that align there (i.e. position != "-")
max_aln_positions = transposed_alignment_df.shape[0] #number of species
percent_aln_positions = transposed_alignment_df.apply(lambda x: len([element for element in x if element != "-"])/max_aln_positions, axis=0) #this returns a series
position_percent_aln_dict = percent_aln_positions.to_dict()

#Select the lowest and highest key which are aligned and in chimeric gene for which value >= 80%
#first and last aligned aminoacid, respectively
chimeric_gene_aln_pos = transposed_alignment_df.loc[my_chimeric_gene,]
chimeric_gene_aligned_pos = list(chimeric_gene_aln_pos[chimeric_gene_aln_pos != "-"].index.values)
valid_positions = [key for key in list(position_percent_aln_dict.keys()) if key in chimeric_gene_aligned_pos]
first_aligned_pos = min([position for position in valid_positions if position_percent_aln_dict[position] >= my_overlap_stringency])
last_aligned_pos = max([position for position in valid_positions if position_percent_aln_dict[position] >= my_overlap_stringency])

#Translate from alignment position to chimeric gene position. Create a dictionary between the two "coordinate" systems (the chimeric gene aln with and without indels)
aln_chimeric_gene_pos_dict = {key: value for value, key in enumerate(valid_positions)}
first_aligned_aa = aln_chimeric_gene_pos_dict[first_aligned_pos] #first aa in the fragment
last_aligned_aa = aln_chimeric_gene_pos_dict[last_aligned_pos] #last aa in the fragment

#Intersect with the aminoacid positions to get the correspondent exon
all_exon_pos_files = glob.glob(my_exon_positions_dir+"/*")
all_exon_pos_df = pd.DataFrame()
for my_file in all_exon_pos_files:
  my_df = pd.read_table(my_file, sep="\t", index_col=False, header=None, names=["geneID", "exon_num","aa_coords", "chr", "genomic_coords", "strand","type"]) #NB: the GeneID col is still in the format ProtID|GeneID
  all_exon_pos_df = pd.concat([all_exon_pos_df, my_df])
all_exon_pos_df["geneID"] = [re.sub(".*\|", "", element) for element in list(all_exon_pos_df["geneID"])] #modify geneID column

#Subset the exon positions dataframe to only the chimeric gene of interest
chimeric_gene_exon_pos = all_exon_pos_df[all_exon_pos_df["geneID"]==my_chimeric_gene]
chimeric_gene_exon_pos_dict = pd.Series(chimeric_gene_exon_pos.exon_num.values, index=chimeric_gene_exon_pos.aa_coords).to_dict() #dictionary with key=aa_pos, value=exon_num
chimeric_gene_exon_pos_ranges = [range(int(element.split("-")[0]), int(element.split("-")[1])) for element in list(chimeric_gene_exon_pos_dict.keys())]
overlapping_exons = [element for element in chimeric_gene_exon_pos_ranges if set(range(first_aligned_aa, last_aligned_aa)).intersection(element)]
first_overlapping_exon = overlapping_exons[0]
last_overlapping_exon = overlapping_exons[len(overlapping_exons)-1]
#Transform ranges back to coordinates
first_overlapping_exon_aa = str(first_overlapping_exon[0])+"-"+str(first_overlapping_exon[len(first_overlapping_exon)-1]+1)
last_overlapping_exon_aa = str(last_overlapping_exon[0])+"-"+str(last_overlapping_exon[len(last_overlapping_exon)-1]+1)
#Intersect with the aminoacid positions to get the correspondent exon numer
first_overlapping_exon_num = chimeric_gene_exon_pos_dict[first_overlapping_exon_aa]
last_overlapping_exon_num = chimeric_gene_exon_pos_dict[last_overlapping_exon_aa]

#Save output file with header:
#1. Chimeric geneID
#2. OG ID
#3. First aligned aa
#4. Last aligned aa
#5. First aligned exon (exon num: aa coords)
#6. Last aligned exon (exon num: aa coords)
final_string = "%s\t%s\t%s\t%s\t%s:%s\t%s:%s" % (my_chimeric_gene, my_orthogroupID, first_aligned_aa, last_aligned_aa, first_overlapping_exon_num, first_overlapping_exon_aa , last_overlapping_exon_num, last_overlapping_exon_aa)
print(final_string)

