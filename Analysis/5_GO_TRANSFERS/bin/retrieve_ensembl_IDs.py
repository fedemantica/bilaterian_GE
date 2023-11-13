#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
import urllib.parse
import urllib.request

parser = argparse.ArgumentParser(description="Script to retrieve the correponding EnsemblID of the provided UniprotKB IDs")
parser.add_argument("--input", "-i", required=True, metavar="input", help="One column files with a list of UniprotKB IDs to translate to EnsemblIDs (no header)")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Two column file containing the UniprotKB ID and the corresponding EnsemblID")
parser.add_argument("--missing", "-m", required=True, metavar="missing", help="One column file containing the UniprotKB ID missing and Ensembl translation")

#Read arguments
args = parser.parse_args()
input_file = args.input
output_file = args.output
missing_file = args.missing

#Generate the query list
input_df = pd.read_table(input_file, sep="\t", index_col=False, header=None, names=["UniprotID"])
my_query_list = list(input_df["UniprotID"])
my_query_string = " ".join(my_query_list) #Generate space-separated string of IDs to be translated. E.g: 'ID1 ID2 ID2'

url = "https://www.uniprot.org/uploadlists/"

params = {
"from": "ACC+ID",
"to": "ENSEMBL_ID",
"format": "tab",
"query": my_query_string
}


original_stdout = sys.stdout

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
   response = f.read()
with open(output_file, "w") as output:
  sys.stdout = output
  print(response.decode('utf-8')) #check what the output looks like
  sys.stdout = original_stdout

#read output file and check which are the missing IDs
translated_ids = list(pd.read_table(output_file, sep="\t", header=0, index_col=False)["From"])
untranslated_ids = pd.Series([element for element in my_query_list if element not in translated_ids]).drop_duplicates()
untranslated_ids.to_csv(missing_file, sep="\t", index=False, header=False, na_rep="NA")
