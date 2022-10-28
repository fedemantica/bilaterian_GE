#!/usr/bin/env python3

import sys
#NB: this has to be run in my Broccoli env to work. Otherwise, decomment the following line (which is ugly)
#sys.path.extend(['/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python36.zip', '/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python3.6', '/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python3.6/lib-dynload', '/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python3.6/site-packages'])

import argparse
from Bio import SeqIO
import pandas as pd
import re
import os
import glob

import subprocess

parser = argparse.ArgumentParser(description="Generate multiple alignments among all the proteins in an orthogroup containing chimeric genes")
parser.add_argument("--input_fastas", "-f", required=True, metavar="input_fastas", help="directory containing all proteome fasta files of species in orthogroups")
parser.add_argument("--input_OG", "-i", required=True, metavar="input_orthogroups", help="orthogroups containing chimeric genes")
parser.add_argument("--chimeric_genes", "-g", required=True, metavar="chimeric_gene", help="comma-separated list of chimeric geneIDs for which to perform the alignment in the relative orthogroup")
parser.add_argument("--mafft", "-m", required=True, metavar="mafft_executable", help="path to mafft execution file")
parser.add_argument("--output", "-o", required=True, metavar="output_dir", help="path to output dir")


#Read arguments
args = parser.parse_args()
my_input_fastas = args.input_fastas
my_input_OG = args.input_OG
my_chimeric_genes = args.chimeric_genes.split(",") #split the original batch
my_mafft = args.mafft
my_output = args.output

print(my_chimeric_genes)
#Header: OG_ID, species, geneID, geneName, chimeric_label
chimeric_OG_df = pd.read_table(str(my_input_OG), sep="\t", header=0, index_col=False)
#Read and join fasta files all together
all_fastas_files = glob.glob(my_input_fastas+"/*") #list all the fasta files in the input directory
all_fastas_entries = []
for my_file in all_fastas_files: #cycle on all the fasta files
  print(my_file)
  my_fasta = list(list(SeqIO.parse(my_file, "fasta")))
  for element in my_fasta: #this I probably don't need if I take the chimeric proteins from the broccoli output 
    element.id = re.sub(".*\|", "", element.id)
    element.description = ""
  all_fastas_entries = all_fastas_entries+my_fasta

for my_chimeric_gene in my_chimeric_genes:
  my_OG_IDs = list(set(list(chimeric_OG_df[chimeric_OG_df["GeneID"]==my_chimeric_gene]["OG_ID"])))
  for my_OG_ID in my_OG_IDs:
    OG_df = chimeric_OG_df[chimeric_OG_df["OG_ID"]==my_OG_ID]
    #remove all the other chimeric genes from the orthogroups
    OG_other_chimeric_genes = [element for element in list(OG_df[OG_df["Chimeric_label"]=="chimeric"]["GeneID"]) if element != my_chimeric_gene]
    OG_genes = [gene for gene in list(OG_df["GeneID"]) if gene not in OG_other_chimeric_genes]
    fastas_entries = [element for element in all_fastas_entries if element.id in OG_genes]
    #save fastas entries to temporary file
    input_temp = my_chimeric_gene+"-"+my_OG_ID+"-input_fasta_tmp.fa"
    with open(input_temp, "w") as input_to_aln:
      SeqIO.write(fastas_entries, input_to_aln, "fasta")
    input_to_aln.close()
    #call shell command to run the mafft alignment
    my_command = "%s --quiet --retree 2 --localpair --maxiterate 1000 input_fasta_tmp.fa > %s/%s-%s-multiple_aln" % (my_mafft, my_output, my_chimeric_gene, my_OG_ID)
    print(my_command) #this is just to know what has been run.
    output_file =  "%s/%s-%s-multiple_aln" % (my_output, my_chimeric_gene, my_OG_ID)
    with open(output_file, "w") as output:
      p = subprocess.Popen([my_mafft, "--quiet", "--retree", "2", "--localpair", "--maxiterate", "1000", input_temp], stdout=output)
      p.wait()
    #remove temporary input
    remove_command = "rm %s" % (input_temp) 
    os.system(remove_command)    
