#!/usr/bin/env python3

import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(description="Script to categorize the chimeric genes based on the overlapping fragments with each orthogroup they belong to")
parser.add_argument("--species", "-s", required=True, metavar="species", help="Species identifier as it appears in the gene orthogroup file")
parser.add_argument("--params_file", "-p", required=True, metavar="species", help="Params file with geneID infos for each species (suffix, length)")
parser.add_argument("--input", "-i", required=True, metavar="input", help="Chimeric orthogroups file")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

#Read arguments
args = parser.parse_args()
my_species = args.species
my_params_file = args.params_file
my_input_file = args.input
my_output_file = args.output

#Read params file
#Header: species, geneID_prefix, geneID_length
params_df = pd.read_table(my_params_file, sep="\t", header=0, index_col=False)
prefix = str(list(params_df[params_df["species"]==my_species]["geneID_prefix"])[0])
length = int(list(params_df[params_df["species"]==my_species]["geneID_length"])[0])

#Read input
#Header: OG_ID, species, geneID, geneName, chimeric_label
input_df = pd.read_table(my_input_file, sep="\t", header=0, index_col=False)
species_chimeric_df = input_df.loc[(input_df["species"]==my_species) & (input_df["chimeric_label"]=="chimeric")]

#Generate the new chimeric geneIDs
chimeric_prefix = prefix+"C"
new_IDs_num = species_chimeric_df.shape[0]*2 #how many new chimeric IDs need to be generated
new_IDs_suffix = [str(element) for element in list(range(1, new_IDs_num+1))]
new_IDs_list = [chimeric_prefix + str("0"*(length-len(chimeric_prefix)-len(element))) + element for element in new_IDs_suffix]
#Pair the chimeric geneIDs
new_IDs_it = iter(new_IDs_list) #NB: the iterator is consumed once called
new_IDs_paired = [element[0]+";"+element[1] for element in list(zip(new_IDs_it, new_IDs_it))]
species_chimeric_df["chimeric_IDs"] = new_IDs_paired
final_df = species_chimeric_df[["geneID", "chimeric_IDs"]] #Add IDs to original df

#Write to output file
final_df.to_csv(my_output_file, sep="\t", index=False, header=True, na_rep="NA")
