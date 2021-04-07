#!/usr/bin/env python3

import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(description="Script to categorize the chimeric genes based on the overlapping fragments with each orthogroup they belong to")
parser.add_argument("--species", "-s", required=True, metavar="species", help="Species identifier as it appears in the gene orthogroup file")
parser.add_argument("--params_file", "-p", required=True, metavar="species", help="Params file with geneID infos for each species (suffix, length)")
parser.add_argument("--input_broken", "-ib", required=True, metavar="input_broken", help="Selected broken genes by species")
parser.add_argument("--input_chimeric", "-ic", required=True, metavar="input_chimeric", help="Chimeric orthogroups file")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

#Read arguments
args = parser.parse_args()
my_species = args.species
my_params_file = args.params_file
my_input_broken = args.input_broken
my_input_chimeric = args.input_chimeric
my_output_file = args.output

#Read params file
#Header: species, geneID_prefix, geneID_length, transcript_suffix, protein_suffix
params_df = pd.read_table(my_params_file, sep="\t", header=0, index_col=False)
prefix = str(list(params_df[params_df["species"]==my_species]["geneID_prefix"])[0])
length = int(list(params_df[params_df["species"]==my_species]["geneID_length"])[0])

########## Broken genes
#Read input
species_broken_df = pd.read_table(my_input_broken, sep="\t", index_col=False, header=None, names=["geneID"])
broken_prefix = prefix+"B"
new_IDs_num = species_broken_df.shape[0] #how many new broken IDs need to be generated
new_IDs_suffix = [str(element) for element in list(range(1, new_IDs_num+1))]
new_IDs_list = [broken_prefix + str("0"*(length-len(broken_prefix)-len(element))) + element for element in new_IDs_suffix]
#Pair the chimeric geneIDs
species_broken_df["new_IDs"] = new_IDs_list #there are only two columns here: broken_geneID, new_IDs
species_broken_df["category"] = "broken"

########## Chimeric genes 
#Read input
#Header: OG_ID, species, geneID, geneName, chimeric_label
chimeric_input_df = pd.read_table(my_input_chimeric, sep="\t", index_col=False, header=0)
chimeric_input_df = chimeric_input_df.loc[(chimeric_input_df["species"]==my_species) & (chimeric_input_df["chimeric_label"]=="chimeric")]
species_chimeric_genes = list(set(list(chimeric_input_df["geneID"])))

#Generate the new chimeric geneIDs
chimeric_prefix = prefix+"C"
new_IDs_num = len(species_chimeric_genes)*2 #how many new chimeric IDs need to be generated
new_IDs_suffix = [str(element) for element in list(range(1, new_IDs_num+1))]
new_IDs_list = [chimeric_prefix + str("0"*(length-len(chimeric_prefix)-len(element))) + element for element in new_IDs_suffix]
#Pair the chimeric geneIDs
new_IDs_it = iter(new_IDs_list) #NB: the iterator is consumed once called
new_IDs_paired = [element[0]+";"+element[1] for element in list(zip(new_IDs_it, new_IDs_it))]

species_chimeric_df = pd.DataFrame({"geneID" : species_chimeric_genes, "new_IDs" : new_IDs_paired, "category" : ["chimeric"]*len(species_chimeric_genes)})

#species_chimeric_df["new_IDs"] = new_IDs_paired
#species_chimeric_df["category"] = "chimeric"
#species_chimeric_df = species_chimeric_df[["geneID", "new_IDs", "category"]] #Add IDs to original df

#concat broken genes info
final_df = pd.concat([species_broken_df, species_chimeric_df])
#Write to output file
final_df.to_csv(my_output_file, sep="\t", index=False, header=True, na_rep="NA")
