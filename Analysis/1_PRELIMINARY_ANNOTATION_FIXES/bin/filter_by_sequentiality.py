#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Script to isolate paralogous genes SEQUENTIALLY annotated on the same chr and on the same strand")
parser.add_argument("--input", "-i", required=True, metavar="input", help="Paralogous genes annotated on the same chr and strand. Format is col1=OG_ID, col2=Species, col3=GeneID")
parser.add_argument("--coord_info", "-c", required=True, metavar="coord_info", help="Genomic start and stop info by gene")
parser.add_argument("--chr_strand_info", "-cs", required=True, metavar="chr_strand_info", help="Chr and strand info by gene")
parser.add_argument("--species", "-s", required=True, metavar="chr_strand_info", help="Species name/identifier as reported in input")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

##### Read arguments #####
args = parser.parse_args()
my_input = args.input
my_coord_info = args.coord_info
my_chr_strand_info = args.chr_strand_info
my_species = args.species
my_output = args.output

##### Main #####
#Read input files
orthogroups_df = pd.read_table(str(my_input), sep="\t", header=None, names=["OG_ID", "Species", "GeneID", "GeneName"])
my_coords_df = pd.read_table(str(my_coord_info), sep="\t", header=None, names=["GeneID", "Chr", "Start", "Stop"])
my_strand_df = pd.read_table(str(my_chr_strand_info), sep="\t", header=None, names=["GeneID", "Chr", "Strand"])
#Create a dictionary with key=GeneID, value=Start
my_start_dict = pd.Series(my_coords_df.Start.values, index=my_coords_df.GeneID).to_dict()
#Create a dictionary with key=GeneID, value=Chr
my_chr_dict = pd.Series(my_coords_df.Chr.values, index=my_coords_df.GeneID).to_dict()
#Create a dictionary with key=GeneID, value=Strand
my_strand_dict = pd.Series(my_strand_df.Strand.values, index=my_strand_df.GeneID).to_dict()

#Initialize final dataframe with the correct columns
final_df = pd.DataFrame(columns=["OG_ID", "Species", "Chr:Strand", "gene1", "gene1_comb", "gene2", "gene2_comb"])
#Group dataframe by OG_ID
my_grouped_df = orthogroups_df.groupby("OG_ID")
for name, group in my_grouped_df:
  genes = list(group["GeneID"])
  #Create a dictionary with key = GeneID, value = (Chr, Strand)
  chr_strand_dict = {element : (my_chr_dict[element], my_strand_dict[element]) for element in genes}
  #Cycle on all possible combinations of paralogous sequentially annotated paralogous genes
  for combination in list(set(list(chr_strand_dict.values()))):
    my_genes = [gene for gene,comb in chr_strand_dict.items() if comb == combination]
    subsetted_start_dict = {key : my_start_dict[key] for key in my_genes}
    ordered_genes = [key for (key, value) in sorted(subsetted_start_dict.items(), key=lambda x: x[1])]
    all_starts = list(my_coords_df.loc[my_coords_df.Chr==my_chr_dict[ordered_genes[0]]]["Start"])
    for i in list(range(len(ordered_genes)-1)):
      first_gene = ordered_genes[i]; first_gene_start=my_start_dict[first_gene]
      second_gene = ordered_genes[i+1]; second_gene_start=my_start_dict[second_gene]
      intermediate_starts = [element for element in all_starts if first_gene_start < element < second_gene_start]
      if len(intermediate_starts) == 0:
        final_entry = pd.DataFrame({"OG_ID" : [name], "Species" : [my_species], "Chr:Strand" : [combination], "gene1" : [first_gene], "gene1_comb" : [chr_strand_dict[first_gene]], "gene2" : [second_gene], "gene2_comb" : [chr_strand_dict[second_gene]]})
        final_df = pd.concat([final_df, final_entry])
#Save to output file
final_df.to_csv(str(my_output), sep="\t", header=False, index=False, na_rep="NA")
