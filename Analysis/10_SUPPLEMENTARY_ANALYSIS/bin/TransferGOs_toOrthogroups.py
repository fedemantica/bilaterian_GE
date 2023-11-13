#!/usr/bin/env python

import argparse  
import pandas as pd
import re

#read arguments
parser = argparse.ArgumentParser(description="Script to transfer the annotation of the query species to all the other species in the same orthogroups whenever number_genes_query/number_genes_withGO <= cutoff")
parser.add_argument("--orthogroups", "-o", required=True, help="Gene orthogroups file with col1=orthogroupID, col2=species, col3=geneID. No header")
parser.add_argument("--GO_file", "-g", required=True, help="GO annotation files with col1=geneID, col2=GOterm, col3=GOterm_description. No header")
parser.add_argument("--species_query", "-s", required=True, help="Species for which GO annotations are provided")
parser.add_argument("--cutoff", "-c", type=int, required=True, help="Cutoff to reduce the GO transfer between onhologs")
parser.add_argument("--output_file", "-out", required=True, help="Output file")

args = parser.parse_args()
orthogroups = args.orthogroups
GO_file = args.GO_file
species_query = args.species_query
cutoff = args.cutoff
output_file = args.output_file

###### Main
#Read files
orthogroups_df = pd.read_table(orthogroups, sep="\t", index_col=False, header=None, names=["ClusterID", "Species", "GeneID"])
orthogroups_df = orthogroups_df[["ClusterID", "Species", "GeneID"]] #subset in case there are more columns in input
GO_df = pd.read_table(GO_file, sep="\t", index_col=False, header=None, names=["GeneID", "GOterm", "GO_description"])
GO_df["GO_annot"] = [re.sub(" ", "_", str(element[0]))+";"+re.sub(" ", "_", str(element[1])) for element in zip(list(GO_df["GOterm"]), list(GO_df["GO_description"]))]

#create a dictionary with key=query_species_gene, value=[GO_term;GO_description, GO_term;GO_description]
grouped_GO_df = GO_df.groupby("GeneID")
collapsed_GO_df = grouped_GO_df.agg({"GO_annot" : lambda x : x.tolist()}) #this is a pandas Series, with index=GeneID
geneID_GO_dict = pd.Series(collapsed_GO_df.GO_annot.values, index=list(collapsed_GO_df.index.values)).to_dict()

#subset the orthogroups dataframe to only query species_genes_genes, translate with corresponding GO list and expand the dataframe
orthogroups_query_df = orthogroups_df.loc[orthogroups_df["Species"]==species_query]
orthogroups_query_df["GO_annot"] = orthogroups_query_df["GeneID"].map(geneID_GO_dict)
orthogroups_query_expanded_df = orthogroups_query_df.explode("GO_annot").reset_index(drop=True)

#group the orthogroups dataframe by clusterID. Select only the GO_term;GO_description entries where support (TOT_query_genes_in_cluster/TOT_query_genes_in_GO) <= cutoff
#create a dictionary with key=clusterID, value=[GO_term;GO_description;support, GO_term;GO_description;support]
grouped_df = orthogroups_query_expanded_df.dropna().groupby("ClusterID") #also remove all columns containing a NA.
filtered_query_species_dict = {}
for name, group in grouped_df:
  filtered_query_species_dict[name] = []
  TOT_genes = len(list(set(list(group["GeneID"]))))
  for GO in list(set(group["GO_annot"])):
    TOT_GO = len(list(set(list(group.loc[group["GO_annot"] == GO]["GeneID"])))) #this we need in case there are duplicated lines in the input file
    support = TOT_genes/TOT_GO
    support_to_print = TOT_GO/TOT_genes
    if support <= cutoff:
      filtered_query_species_dict[name].append(GO+";"+str(support_to_print))

#Add column with [GO_term;GO_description;support] to the original orthogroups dataframe
orthogroups_df["GO_annot"] = orthogroups_df["ClusterID"].map(filtered_query_species_dict)
#Expand the original orthogroups dataframe dataframe
orthogroups_GO_df = orthogroups_df.explode("GO_annot").reset_index(drop=True)
#Separate the last entry into three different fields and remove original one
orthogroups_GO_df[["GO_term","GO_description","Support"]] = orthogroups_GO_df["GO_annot"].str.split(";", expand=True)
orthogroups_GO_df = orthogroups_GO_df.drop(columns="GO_annot")
#Save to file
orthogroups_GO_df.to_csv(output_file, sep="\t", header=True, index=False, na_rep="NA")
