#Rscript to simply compute the log@ of the input table
Args = commandArgs(trailingOnly=TRUE);
TPM_file = (Args[1])
output_file = (Args[2])

#Upload libraries
#library("tidyverse")

#Read arguments in
input_df = read.delim(TPM_file, header=TRUE, sep="\t", row=1)

#Compute log2
log2_df = apply(input_df, c(1,2), function(x) log2(x+1))
 
#Save to output file
write.table(log2_df, output_file, quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)