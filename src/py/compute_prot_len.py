#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import re

parser = argparse.ArgumentParser(description="Script to compute protein lenght starting from a protein fasta file")
parser.add_argument("--input", "-i", required=True, metavar="input_file", help="Protein fasta file")
parser.add_argument("--output", "-o", required=True, metavar="output_file", help="Two columns file contaning with col1=ProteinID, col2=ProteinLength")

#Read arguments
args = parser.parse_args()
input_file = args.input
output_file = args.output

#Open output filehandle
my_output = open(str(output_file), "w")
#Read fasta input
records = list(SeqIO.parse(str(input_file), "fasta"))
for element in records:
  my_id = re.sub(".*\\|", "", element.id)
  my_len = len(str(element.seq)) #compute length of the protein
  my_output.write("%s\t%d\n" % (my_id, my_len))
