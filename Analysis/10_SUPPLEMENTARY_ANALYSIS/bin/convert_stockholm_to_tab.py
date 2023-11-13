#!/usr/bin/python
import sys
import re

if len(sys.argv) != 3:
    print("Usage: python convert_sto_to_tsv.py input.sto output.tsv")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

try:
    with open(input_file, "rb") as infile, open(output_file, "w") as outfile:
        seq_id = ""
        clan_id = ""
        mb_values = []

        for line in infile:
            line = line.strip()

            if line.startswith(b"#=GF ID"):
                seq_id = line.split()[2]
            elif line.startswith(b"#=GF AC"):
                clan_id = line.split()[2]
            elif line.startswith(b"#=GF MB"):
                mb_value = line.split()[2]
                mb_values.append(mb_value)

            elif line.startswith(b"//") and seq_id and mb_values:
                # Output the required information in the desired format
                for mb_value in mb_values:
                    final_string = str(mb_value)+"\t"+str(clan_id)+"\t"+str(seq_id)+"\n"
                    formatted_final_string = re.sub(";", "", re.sub("'", "", re.sub("b'", "", final_string)))
                    outfile.write(formatted_final_string)

                # Reset variables for the next entry
                seq_id = ""
                clan_id = ""
                mb_values = []

except FileNotFoundError:
    print("Input file not found:", input_file)
    sys.exit(1)
