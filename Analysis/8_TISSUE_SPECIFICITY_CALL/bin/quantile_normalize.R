#General Rscript which runs normalizeBetweenArrays from Limma on the input table
.libPaths(c("/nfs/users/mirimia/fmantica/software/anaconda3/envs/R3.6_env/lib", .libPaths()))

library("limma")

#allow external arguments. Set input and output.
Args <- commandArgs(trailingOnly=TRUE);
input_file = (Args[1]);
output_file = (Args[2]);

input_table = read.table(input_file, header=TRUE, row=1)
normalized_table = normalizeBetweenArrays(input_table, method="quantile")
write.table(normalized_table, file=output_file, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
