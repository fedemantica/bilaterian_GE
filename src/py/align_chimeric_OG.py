#!/usr/bin/env python3

import sys
#NB: this has to be run in my Broccoli env to work. Otherwise, decomment the following line (which is ugly)
#sys.path.extend(['/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python36.zip', '/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python3.6', '/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python3.6/lib-dynload', '/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python3.6/site-packages'])

import argparse
from Bio import SeqIO
import pandas as pd
import re
import glob
import os

parser = argparse.ArgumentParser(description="Generate multiple alignments among all the proteins in an orthogroup containing chimeric genes")
parser.add_argument("--input_fastas", "-f", required=True, metavar="input_fastas", help="directory containing all proteome fasta files of species in orthogroups")
parser.add_argument("--input_OG", "-i", required=True, metavar="input_orthogroups", help="orthogroups containing chimeric genes")
parser.add_argument("--chimeric_gene", "-g", required=True, metavar="chimeric_gene", help="geneID of a chimeric gene as reported in the fasta files")
parser.add_argument("--mafft", "-m", required=True, metavar="mafft_executable", help="path to mafft execution file")
parser.add_argument("--output", "-o", required=True, metavar="output_file", help="path to output file")


#Read arguments
args = parser.parse_args()
my_input_fastas = args.input_fastas
my_input_OG = args.input_OG
my_chimeric_gene = args.chimeric_gene
my_mafft = args.mafft
my_output = args.output

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

my_OG_ID = list(chimeric_OG_df[chimeric_OG_df["geneID"]==my_chimeric_gene]["OG_ID"])[0]
OG_df = chimeric_OG_df[chimeric_OG_df["OG_ID"]==my_OG_ID]
#I still need to remove all the other chimeric genes from the orthogroup.
OG_genes = list(OG_df["geneID"])
fastas_entries = [element for element in all_fastas_entries if element.id in OG_genes]
#save fastas entries to temporary file
with open('input_fasta_tmp.fa', 'w') as input_to_aln:
  SeqIO.write(fastas_entries, input_to_aln, "fasta")
#call shell command to run the mafft alignment
my_command = "%s --quiet --retree 2 --localpair --maxiterate 1000 input_fasta_tmp.fa > %s" % (my_mafft, my_output)
print(my_command)
os.system(my_command)
