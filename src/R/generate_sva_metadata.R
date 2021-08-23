#Rscript to be run from a snakemake
#It generates the metadata tables needed to run the sva correction
#Probably to be improved when I modify the function to compute the inter-tissue correlation.

Args = commandArgs(trailingOnly=TRUE);
metadata_input = (Args[1]) 
metadata_output = (Args[2]) 

#metadata_input = "/Users/federica/mnt/projects/bilaterian_GE/data/samples_metadata/Hs2_samples_info.tab"

metadata_df = read.delim(metadata_input, header=TRUE)
metadata_df = metadata_df[,c("Species", "Tissue", "Sample", "Read_number")]
colnames(metadata_df) = c("scientific_name", "tissue", "run", "num_read_unfiltered")

#change type of num_read_unfiltered column
metadata_df$num_read_unfiltered = as.numeric(gsub(",", "", as.vector(metadata_df$num_read_unfiltered)))

#add new mock required fields
metadata_df$exclusion = rep("no", nrow(metadata_df))
metadata_df$bioproject = seq(1, nrow(metadata_df))
metadata_df$num_read_fastp = metadata_df$num_read_unfiltered
metadata_df$num_read_masked = metadata_df$num_read_unfiltered
metadata_df$mapping_rate_masked = rep(0.8, nrow(metadata_df))
metadata_df$instrument = rep("not_provided", nrow(metadata_df))

#change data type
metadata_df$tissue = as.vector(metadata_df$tissue)

#saving to output
write.table(metadata_df, metadata_output, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)