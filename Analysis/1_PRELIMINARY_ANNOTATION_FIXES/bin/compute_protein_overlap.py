#!/usr/bin/env python3

import sys
#NB: this has to be run in my Broccoli env to work. Otherwise, decomment the following line (which is ugly)
sys.path.extend(['/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python36.zip', '/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python3.6', '/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python3.6/lib-dynload', '/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python3.6/site-packages'])

import argparse
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Compute the percentage of overlap (i.e. not gapped positions) between two aligned proteins")
parser.add_argument("--alignment", "-i", required=True, metavar="alignment_file", help="fasta file with the alignes proteins")
#parser.add_argument("--output", "-o", required=True, metavar="output_file")

#overlap is the sequence of the protein which is matched by the other one. It is defined as the difference between the alignment length and the number of gaps in the OTHER protein alignment, divided by the alignment length".
def compute_overlap(records):
  gene1 = records[0].id
  gene2 = records[1].id
  if len(str(records[0].seq)) == len(str(records[1].seq)):
    aln_len = len(str(records[0].seq))
    #I can't stress enough the fact that I am looking at the gaps in the OTHER protein
    overlap_gene1 = (aln_len - len([element for element in str(records[1].seq) if element=="-"]))/aln_len
    overlap_gene2 = (aln_len - len([element for element in str(records[0].seq) if element=="-"]))/aln_len
    print("%s\t%s\t%f" % (gene1, gene2, overlap_gene1)) #I am breaking this so that it is more readable
    print("%s\t%s\t%f" % (gene2, gene1, overlap_gene2))
  else:
    raise Exception("The input is not a proper alignment: different lenghts for the two proteins")

def main():
  args = parser.parse_args() #read arguments
  #my_output = open(args.output, 'w') #open output filehandle
  records = list(list(SeqIO.parse(args.alignment, "fasta")))
  compute_overlap(records)

main()
