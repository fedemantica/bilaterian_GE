#!/usr/bin/env python

import argparse  
import pandas as pd
import re

#read arguments
parser = argparse.ArgumentParser(description="Script to transfer the annotation of the query species to all the other species in the same orthogroups whenever number_genes_query/number_genes_withGO <= cutoff")
parser.add_argument("--orthogroups", "-o", required=True, help="Gene orthogroups file with col1=orthogroupID, col2=species, col3=geneID. No header")
parser.add_argument("--clade_specific_GO_file", "-cs", required=True, help="Clade specific GO annotation file with col1=GOterm, col2=GOterm_description. No header")
parser.add_argument("--GO_file", "-g", nargs="+", required=True, help="GO annotation files with col1=geneID, col2=GOterm, col3=GOterm_description. No header")
parser.add_argument("--cutoff", "-c", type=int, required=True, help="Cutoff to reduce the GO transfer between onhologs")
parser.add_argument("--output_file", "-out", required=True, help="Output file")

args = parser.parse_args()
orthogroups = args.orthogroups
clade_specific_GO_file = args.clade_specific_GO_file
GO_files_list = args.GO_file #This argument is a list
cutoff = args.cutoff
output_file = args.output_file

###### Main
#Read files
orthogroups_df = pd.read_table(orthogroups, sep="\t", index_col=False, header=None, names=["ClusterID", "Species", "GeneID"])
orthogroups_df = orthogroups_df[["ClusterID", "Species", "GeneID"]] #subset in case there are more columns in input
clade_specific_GO_df = pd.read_table(clade_specific_GO_file, sep="\t", index_col=False, header=None, names=["GOterm", "GO_description"])

GO_df = pd.DataFrame()
for GO_file in GO_files_list:
  species_GO_df = pd.read_table(GO_file, sep="\t", index_col=False, header=None, names=["GeneID", "GOterm", "GO_description"])
  species_GO_df["GO_annot"] = [re.sub(" ", "_", str(element[0]))+";"+re.sub(" ", "_", str(element[1])) for element in zip(list(species_GO_df["GOterm"]), list(species_GO_df["GO_description"]))]
  GO_df = pd.concat([GO_df, species_GO_df])

#Filter for only the GOs that are clade specific
GO_df = GO_df.loc[GO_df["GOterm"].isin(list(clade_specific_GO_df["GOterm"]))]

#create a dictionary with key=query_species_gene, value=[GO_term;GO_description, GO_term;GO_description]
grouped_GO_df = GO_df.groupby("GeneID")
collapsed_GO_df = grouped_GO_df.agg({"GO_annot" : lambda x : x.tolist()}) #this is a pandas Series, with index=GeneID
geneID_GO_dict = pd.Series(collapsed_GO_df.GO_annot.values, index=list(collapsed_GO_df.index.values)).to_dict()

#Translate the geneID with the corresponding GO list and expand the dataframe
orthogroups_query_df = orthogroups_df.copy()
orthogroups_query_df["GO_annot"] = orthogroups_query_df["GeneID"].map(geneID_GO_dict)
orthogroups_query_expanded_df = orthogroups_query_df.explode("GO_annot").reset_index(drop=True)

#group the orthogroups dataframe by clusterID. Select only the GO_term;GO_description entries where support (TOT_query_genes_in_cluster/TOT_query_genes_in_GO) <= cutoff FOR ALL SPECIES
#create a dictionary with key=clusterID, value=[GO_term;GO_description;support, GO_term;GO_description;support]
grouped_df = orthogroups_query_expanded_df.dropna().groupby(["ClusterID", "GO_annot"]) #also remove all columns containing a NA.
filtered_query_species_dict = {}
for name, group in grouped_df:
  filtered_query_species_dict[name[0]] = []

for name, group in grouped_df:
  ClusterID = name[0]
  GO = name[1] #Get the GO
  all_supports = []
  all_supports_to_print = []
  #filtered_query_species_dict[ClusterID] = []
  subgroup_df = group.groupby("Species")
  for subname, subgroup in subgroup_df:
    TOT_genes = len(list(set(list(orthogroups_df.loc[(orthogroups_df["ClusterID"]==ClusterID) & (orthogroups_df["Species"]==subname)])))) #Total number of genes from a species in that orthogroup
    TOT_GO = len(list(set(list(subgroup["GeneID"])))) #Total number of genes of a species in the tested GO
    support = TOT_genes/TOT_GO
    support_to_print = TOT_GO/TOT_genes
    all_supports_to_print = all_supports_to_print + [support_to_print]
    all_supports = all_supports + [support]
  if all(element <= cutoff for element in all_supports):
    filtered_query_species_dict[ClusterID].append(GO+";"+str(max(all_supports_to_print)))

#Add column with [GO_term;GO_description;support] to the original orthogroups dataframe
orthogroups_df["GO_annot"] = orthogroups_df["ClusterID"].map(filtered_query_species_dict)
#Expand the original orthogroups dataframe dataframe
orthogroups_GO_df = orthogroups_df.explode("GO_annot").reset_index(drop=True)
#Separate the last entry into three different fields and remove original one
orthogroups_GO_df[["GO_term","GO_description","Support"]] = orthogroups_GO_df["GO_annot"].str.split(";", expand=True)
orthogroups_GO_df = orthogroups_GO_df.drop(columns="GO_annot")
#Save to file
orthogroups_GO_df.to_csv(output_file, sep="\t", header=True, index=False, na_rep="NA")
