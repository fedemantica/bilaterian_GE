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
import time
import numpy as np

parser = argparse.ArgumentParser(description="Script to recluster the gene orthogroups based on the orthopairs connections between a subset of species")
parser.add_argument("--orthogroups", "-og", required=True, metavar="orthogroups", help="Corrected gene orthogroups, with col1=orthogroupID, col2=species, col3=geneID")
parser.add_argument("--score_matrix", "-a", required=True, metavar="score_matrix", help="Score matrix in the long format: col1=AA1, col2=AA2, col3=Score")
parser.add_argument("--OG_ID", "-id", required=True, metavar="OG_ID", help="OrthogroupID of the orthogroup of interest")
parser.add_argument("--input_fastas", "-f", required=True, metavar="input_fastas", help="Directory containing all proteome fasta files of species in orthogroups")
parser.add_argument("--mafft", "-m", required=True, metavar="mafft_executable", help="path to mafft execution file")
parser.add_argument("--gap_penalty", "-g", required=True, metavar="gap_penalty", help="penalty for gap opening")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output directory")

###### Read arguments
args = parser.parse_args()
orthogroups_file = args.orthogroups
score_file = args.score_matrix
my_orthogroup = args.OG_ID
input_fastas_dir = args.input_fastas
my_mafft = args.mafft
gap_penalty = args.gap_penalty
output_dir = str(args.output)

##################################
###### DEFINE FUNCTIONS ##########
##################################

def compute_sim_score(aln, AA_pair_score_dict):
  first_gene = aln[0].seq
  second_gene = aln[1].seq
  final_score = 0
  for position in list(range(0, len(first_gene))):
    AA_pair = str(first_gene[position]) + str(second_gene[position])
    score = AA_pair_score_dict[AA_pair]
    final_score = final_score + score
  first_gene_score = final_score/len(re.sub("-", "", str(first_gene)))
  second_gene_score = final_score/len(re.sub("-", "", str(second_gene)))
  return([first_gene_score, second_gene_score])

##################################
###### READ INPUTS ###############
##################################
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

#Read dataframe with scores for the AA pairs
score_df = pd.read_table(str(score_file), sep="\t", index_col=False, header=None, names=["First_AA", "Second_AA", "Score"])
score_df["AA_pair"] = score_df["First_AA"]+score_df["Second_AA"]
#Generate dictionary with key=AA_pair and value=Score
AA_pair_score_dict = pd.Series(score_df.Score.values, index=score_df.AA_pair).to_dict()
#Add entries to dictionary
all_indel_pairs = ["-"+element for element in list(set(score_df["First_AA"]))] + [element+"-" for element in list(set(score_df["First_AA"]))]
for indel_pair in all_indel_pairs:
  AA_pair_score_dict[indel_pair] = 0 #count zero if there is a mismatch

##################################
###### MAIN ######################
##################################
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
  #Change ids of fastas entries to geneID1_geneID2
  fastas_entries[0].id = geneID1+"_"+geneID2
  fastas_entries[1].id = geneID2+"_"+geneID1
  #Save fastas entries to temporary file
  input_temp = output_dir+"/"+geneID1+"-"+geneID2+"-input_fasta_tmp.fa"
  with open(input_temp, "w") as input_to_aln:
    SeqIO.write(fastas_entries, input_to_aln, "fasta")
  input_to_aln.close()
  #Write shell command to print to standard output
  my_command = "%s --quiet --retree 2 --localpair --maxiterate 1000 input_fasta_tmp.fa > %s/%s-%s-%s-pairwise_aln" % (my_mafft, output_dir, geneID1, geneID2, my_orthogroup)
  print(my_command) #this is for debugging
  output_aln_file = "%s/%s-%s-%s-pairwise_aln" % (output_dir, geneID1, geneID2, my_orthogroup)
  #write to output file
  with open(output_aln_file, "w") as output_aln:
    p = subprocess.Popen([my_mafft, "--quiet", "--retree", "2", "--localpair", "--maxiterate", "1000", input_temp], stdout=output_aln)
 
  #This is necessary because the I can't copy the SeqIO object
  fastas_entries[0].id = geneID1
  fastas_entries[1].id = geneID2

  #compute the sequence similarity
  #There are some pairs of AA for which we just count 1.
  time.sleep(0.3) #Delay to have time to write to output file
  aln_output = list(list(SeqIO.parse(output_aln_file, "fasta"))) #If not enough, delay more.
  if len(aln_output) == 0:
    time.sleep(1)
    aln_output = list(list(SeqIO.parse(output_aln_file, "fasta")))
    if len(aln_output) == 0:
      time.sleep(1.5)
      aln_output = list(list(SeqIO.parse(output_aln_file, "fasta")))

  sim_score = compute_sim_score(aln_output, AA_pair_score_dict) #This is a list containing the sim_score of geneID1 and the sim_score of geneID2
  output_sim_score_file = "%s/%s-%s-%s-sim_scores" % (output_dir, geneID1, geneID2, my_orthogroup)
  #Print the similarity score in a reciprocal way.
  with open(output_sim_score_file, "w") as output_sim_score:
    output_sim_score.write("%s\t%s\t%s\t%f\n" % (my_orthogroup, geneID1, geneID2, sim_score[0]))
    output_sim_score.write("%s\t%s\t%s\t%f\n" % (my_orthogroup, geneID2, geneID1, sim_score[1]))
  output_sim_score.close()

#Cat all files from the single alignments into a final one
cat_aln_command = "cat %s/*-%s-pairwise_aln > %s/%s-all_pairwise_aln" % (output_dir, my_orthogroup, output_dir, my_orthogroup)
cat_score_command = "cat %s/*-%s-sim_scores > %s/%s-all_sim_scores" % (output_dir, my_orthogroup, output_dir, my_orthogroup)
os.system(cat_aln_command)
os.system(cat_score_command)
#Generate commands to remove the input and not-concatenated files
single_pairwise_rm_command = "rm %s/*-%s-pairwise_aln" % (output_dir, my_orthogroup)
single_sim_scores_rm_command = "rm %s/*-%s-sim_scores" % (output_dir, my_orthogroup)
input_rm_command = "rm %s/*-input_fasta_tmp.fa" % (output_dir)
#Execute remove commands
os.system(single_pairwise_rm_command)
os.system(single_sim_scores_rm_command)
os.system(input_rm_command)

#Re-read input and compute the average sequence similarity for each gene.
#The output should have: OG_ID, Species, GeneID, average_similarity
all_sim_scores_file = "%s/%s-all_sim_scores" % (output_dir, my_orthogroup)
all_sim_scores_df = pd.read_table(all_sim_scores_file, sep="\t", index_col=False, header=None, names=["OG_ID", "GeneID1", "GeneID2", "Sim_score"])
all_sim_scores_grouped_df = all_sim_scores_df.groupby("GeneID1")
average_sim_df = pd.DataFrame()
for gene, group in all_sim_scores_grouped_df:
  average_sim_score = np.mean(list(group["Sim_score"]))
  species = geneID_species_dict[gene]
  average_sim_df = pd.concat([average_sim_df, pd.DataFrame({"OG_ID" : [my_orthogroup], "Species" : [species], "GeneID" : [gene], "Avg_sim_score" : [average_sim_score]})])
average_sim_file =  "%s/%s-average_sim_scores" % (output_dir, my_orthogroup)
average_sim_df.to_csv(average_sim_file, sep="\t", index=False, header=False, na_rep="NA")
