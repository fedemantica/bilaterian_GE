#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
import re
import os
import glob
import itertools
from Bio import SeqIO
import subprocess

parser = argparse.ArgumentParser(description="Script to recluster the gene orthogroups based on the orthopairs connections between a subset of species")
parser.add_argument("--orthogroups", "-og", required=True, metavar="orthogroups", help="Corrected gene orthogroups, with col1=orthogroupID, col2=species, col3=geneID")
parser.add_argument("--OG_ID", "-id", required=True, metavar="orthogroups", help="OrthogroupID of the orthogroup of interest")
parser.add_argument("--input_fastas", "-f", required=True, metavar="input_fastas", help="Directory containing all proteome fasta files of species in orthogroups")
parser.add_argument("--mafft", "-m", required=True, metavar="mafft_executable", help="path to mafft execution file")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output directory")

###### Read arguments
args = parser.parse_args()
orthogroups_file = args.orthogroups
my_orthogroup = args.OG_ID
input_fastas_dir = args.input_fastas
my_mafft = args.mafft
output_dir = str(args.output)

##################################
###### MAIN ######################
##################################
### Read inputs
#Orthogroups
orthogroups_df = pd.read_table(orthogroups_file, sep="\t", index_col=False, header=None, names=["OG_ID", "Species", "GeneID"])
#Generate geneID-species dict
geneID_species_dict = pd.Series(orthogroups_df.Species.values, index=orthogroups_df.GeneID).to_dict()

#Join fasta files all together
all_fastas_files = glob.glob(input_fastas_dir+"/*") #list all the fasta files in the input directory
all_fastas_entries = []
for my_file in all_fastas_files: #cycle on all the fasta files
  print(my_file)
  my_fasta = list(list(SeqIO.parse(my_file, "fasta")))
  for element in my_fasta: 
    element.id = re.sub(".*\|", "", element.id)
    element.description = ""
  all_fastas_entries = all_fastas_entries+my_fasta


#Select the orthogroup of interest
filtered_orthogroup_df = orthogroups_df.loc[orthogroups_df["OG_ID"]==my_orthogroup]

#Generate all possible pair of genes between genes belonging
all_genes = list(filtered_orthogroup_df["GeneID"])
all_gene_pairs = list(itertools.combinations(all_genes, 2))
#Filter out pairs where genes come from the same species
filt_gene_pairs = [gene for gene in all_gene_pairs if geneID_species_dict[gene[0]] != geneID_species_dict[gene[1]]]
#Transform tuple to string: [("geneID1", "geneID2")] to ["geneID1_geneID2"]
#filt_gene_pairs = [gene[0]+"_"+gene[1] for gene in filt_gene_pairs]
#Transform tuple to a list
filt_gene_pairs = [[gene[0], gene[1]] for gene in filt_gene_pairs]

#Align each pair of genes and save alignment to output
for gene_pair in filt_gene_pairs:
  geneID1 = gene_pair[0]
  geneID2 = gene_pair[1]
  #This is to maintain the order of the IDs
  fastas_entries_1 = [element for element in all_fastas_entries if element.id==geneID1][0]
  fastas_entries_2 = [element for element in all_fastas_entries if element.id==geneID2][0]
  fastas_entries = [fastas_entries_1, fastas_entries_2]
  print(len(fastas_entries))
  print(type(fastas_entries))
  print(fastas_entries)
  #geneID1_original = str(fastas_entries[0].id)
  #geneID2_original = str(fastas_entries[1].id)
  #Change ids of fastas entries to geneID1_geneID2
  #geneID1_new_ID = fastas_entries[0].id+"_"+fastas_entries[1].id
  #geneID2_new_ID = fastas_entries[1].id+"_"+fastas_entries[0].id
  fastas_entries[0].id = geneID1+"_"+geneID2
  fastas_entries[1].id = geneID2+"_"+geneID1
  print(fastas_entries)
  #Save fastas entries to temporary file
  input_temp = output_dir+"/"+geneID1+"-"+geneID2+"-input_fasta_tmp.fa"
  with open(input_temp, "w") as input_to_aln:
    SeqIO.write(fastas_entries, input_to_aln, "fasta")
  input_to_aln.close()
  #Write shell command to print to standard output
  my_command = "%s --quiet --retree 2 --localpair --maxiterate 1000 input_fasta_tmp.fa > %s/%s-%s-%s-pairwise_aln" % (my_mafft, output_dir, geneID1, geneID2, my_orthogroup)
  print(my_command) #this is for debugging
  output_file = "%s/%s-%s-%s-pairwise_aln" % (output_dir, geneID1, geneID2, my_orthogroup)
  #write to output file
  with open(output_file, "w") as output:
    p = subprocess.Popen([my_mafft, "--quiet", "--retree", "2", "--localpair", "--maxiterate", "1000", input_temp], stdout=output)
  #This is necessary because the I can't copy the SeqIO object
  fastas_entries[0].id = geneID1
  fastas_entries[1].id = geneID2
  #remove temporary input
  #remove_command = "rm %s" % (input_temp)
  #os.system(remove_command)

#Cat all files from the single alignments into a final one
cat_command = "cat %s/*-%s-pairwise_aln > %s/%s-all_pairwise_aln" % (output_dir, my_orthogroup, output_dir, my_orthogroup)
os.system(cat_command)
#Remove the single alignment files
single_pairwise_remove_command = "rm %s/*-%s-pairwise_aln" % (output_dir, my_orthogroup)
input_remove_command = "rm %s/*-input_fasta_tmp.fa" % (output_dir)
os.system(single_pairwise_remove_command)
os.system(input_remove_command)
