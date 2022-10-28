##### Rscript to compute transcript and gene-wise TPMs and RPKMs starting from normalized counts (with DESeq2)

##### Update .libPaths()
#.libPaths(c("/Users/federica/mnt/software/R/3.5", "/Users/federica/mnt/software/R", .libPaths()))

library("DESeq2")
library("hashmap")
library("tidyverse")
library("tximportData")
library("tximport")

### Read arguments
Args = commandArgs(trailingOnly=TRUE);
species = (Args[1])
tx2gene_file = (Args[2])
transcript_counts_file = (Args[3])
metadata_file = (Args[4])
samples_folder = (Args[5])
output_folder = (Args[6])

print(species)
print(tx2gene_file)
print(transcript_counts_file)
print(metadata_file)
print(samples_folder)
print(output_folder)

#For debugging
# home = "/Users/federica/mnt/"
# species = "Mdo"
# tx2gene_file = paste0(home, "projects/bilaterian_GE/data/preprocessing_all/kallisto_out/Mdo/transcript_gene_dict.tab")
# transcript_counts_file = paste0(home, "projects/bilaterian_GE/data/preprocessing_all/kallisto_out/Mdo/all_samples_transcript_counts.tab")
# metadata_file = paste0(home, "projects/bilaterian_GE/data/samples_metadata/Mdo_samples_info.tab")
# samples_folder = paste0(home, "projects/bilaterian_GE/data/preprocessing_all/kallisto_out/Mdo/")
# output_folder = paste0(home, "projects/mixed_files/TPM_expr_table/Mdo_test/")

## Upload transcript-count table
transcript_counts_table = read.delim(transcript_counts_file, header=TRUE, row=1)
#Round counts
transcript_counts_table = apply(transcript_counts_table, c(1,2), round)
ordered_samples = colnames(transcript_counts_table)

## Upload metadata table
metadata_df = unique(read.delim(metadata_file, header=TRUE)[,1:4])
colnames(metadata_df) = c("species", "tissues", "metasample", "sample")
metadata_df$species = rep(species, nrow(metadata_df))
rownames(metadata_df) = as.vector(metadata_df$sample)
metadata_df$metasample = NULL
metadata_df$sample = NULL
#order as in the transcript counts table
metadata_df = metadata_df[ordered_samples,]

## Generate list of all samples original kallisto out
all_files = list.files(path=samples_folder, pattern="abundance.tsv", recursive = TRUE, full.names=TRUE)

### Normalize counts with DESeq2
#Create DESeq2 object. No design
dds <- DESeqDataSetFromMatrix(countData = as.matrix(transcript_counts_table), colData = metadata_df, design = ~1)
dds <- estimateSizeFactors(dds)
normalized_counts_df <- as.data.frame(counts(dds, normalized=TRUE))

#Cycle on each sample
ordered_genes = rownames(normalized_counts_df) 
ALL_tpm_transcript_table = data.frame(matrix(nrow = nrow(normalized_counts_df), 
                                         ncol=0, 
                                         dimnames = list(ordered_genes, NULL)))
ALL_rpkm_transcript_table = ALL_tpm_transcript_table

for (sample in colnames(normalized_counts_df)) {
  ## Upload original abundance.tsv file
  #isolate original file from all_files listing
  sample_file = all_files[grep(paste0(sample, "/abundance.tsv"), all_files)]
  sample_original_df = read.delim(sample_file, header=TRUE, row=1)
  
  ## TPMs
  #### Take the first transcript with counts != 0
  selected_transcript_entry = subset(sample_original_df, est_counts > 0)[1,]
  ### Compute scaling factor (SF) from original data
  tpm_SF = selected_transcript_entry$est_counts / (selected_transcript_entry$eff_length * selected_transcript_entry$tpm)
  ### Compute transcript-level TPMs with normalized counts
  sample_norm_counts_df = normalized_counts_df[,sample, drop=FALSE] #select sample from the whole normalized counts table
  eff_length_df = sample_original_df[,"eff_length", drop=FALSE] #select the effective length from the original kallisto output
  joint_df = merge(sample_norm_counts_df, eff_length_df, by="row.names") #merge the tables with the norm counts and effective length
  colnames(joint_df) = c("transcriptID", "norm_counts", "eff_length") #update colnames
  joint_df$tpm_norm_counts = (joint_df$norm_counts / joint_df$eff_length) / tpm_SF #actually compute the tpms
  rownames(joint_df) = joint_df$transcriptID #update rownames with transcriptID
  tpm_transcript_table = joint_df[ordered_genes, "tpm_norm_counts", drop=FALSE] #isolate only the norm_tpm info
  colnames(tpm_transcript_table) = sample #add sample as norm_tpm name
  ALL_tpm_transcript_table = cbind(ALL_tpm_transcript_table, tpm_transcript_table) #join the sample info to the final table.
  
  ## RPMKs
  ### Compute scaling factor (SF): sum of all counts / 10^6
  tot_counts = sum(sample_original_df$est_counts) #get the sum of all counts
  rpkm_SF = tot_counts/10^6 #compute the RPKM scaling factor
  ### Compute transcript-level RPKMs with normalized counts
  joint_df$rpkm_norm_counts = (joint_df$norm_counts / rpkm_SF) / joint_df$eff_length #actually compute the rpkms
  rpkm_transcript_table = joint_df[ordered_genes, "rpkm_norm_counts", drop=FALSE] #isolate only the norm_rpkm info
  colnames(rpkm_transcript_table) = sample #add sample as norm_rpkm name
  ALL_rpkm_transcript_table = cbind(ALL_rpkm_transcript_table, rpkm_transcript_table)
}

### Generate transcriptID : geneID dictionary
tx2gene_table = read.delim(tx2gene_file, header=TRUE)
tx2gene_dict = hashmap(as.vector(tx2gene_table$transcriptID), as.vector(tx2gene_table$geneID))
### Add geneID as column to transcript quantification tables (tpm or crpkms)
ALL_tpm_gene_table = ALL_tpm_transcript_table
ALL_rpkm_gene_table = ALL_rpkm_transcript_table
ALL_tpm_gene_table$geneID = tx2gene_dict$find(rownames(ALL_tpm_transcript_table))
ALL_rpkm_gene_table$geneID = tx2gene_dict$find(rownames(ALL_rpkm_transcript_table))

### Compute gene-level TPMs (sum across transcript-level quantification)
#I am not completely sure how this function can work, but it really does.
ALL_tpm_gene_table_final = aggregate(ALL_tpm_gene_table[,colnames(ALL_tpm_gene_table)[colnames(ALL_tpm_gene_table) != "geneID"]], 
                                     by=list(ALL_tpm_gene_table$geneID), 
                                     sum)
#Format dataframe
colnames(ALL_tpm_gene_table_final)[1] = "geneID"
rownames(ALL_tpm_gene_table_final) = ALL_tpm_gene_table_final$geneID
ALL_tpm_gene_table_final$geneID = NULL

### Compute gene-level TPMs (sum across transcript-level quantification)
ALL_rpkm_gene_table_final = aggregate(ALL_rpkm_gene_table[,colnames(ALL_rpkm_gene_table)[colnames(ALL_rpkm_gene_table) != "geneID"]], 
                                     by=list(ALL_rpkm_gene_table$geneID), 
                                     sum)
#Format dataframe
colnames(ALL_rpkm_gene_table_final)[1] = "geneID"
rownames(ALL_rpkm_gene_table_final) = ALL_rpkm_gene_table_final$geneID
ALL_rpkm_gene_table_final$geneID = NULL

### Compute the log2
ALL_tpm_gene_table_log2 = apply(ALL_tpm_gene_table_final, c(1,2), function(x) log2(x+1)) 
ALL_rpkm_gene_table_log2 = apply(ALL_rpkm_gene_table_final, c(1,2), function(x) log2(x+1))

### Save all tables to output file
#norm transcript counts
write.table(normalized_counts_df, paste0(output_folder, "all_samples_transcript-normcounts.tab"), quote=FALSE, sep="\t", col.names=NA)
#norm transcript tpms
write.table(ALL_tpm_transcript_table, paste0(output_folder, "all_samples_transcript_TPMs-normcounts.tab"), quote=FALSE, sep="\t", col.names=NA)
#norm transcript rpkms
write.table(ALL_rpkm_transcript_table, paste0(output_folder, "all_samples_transcript_RPKMs-normcounts.tab"), quote=FALSE, sep="\t", col.names=NA)
#norm gene tpms
write.table(ALL_tpm_gene_table_final, paste0(output_folder, "all_samples_gene_TPMs-normcounts.tab"), quote=FALSE, sep="\t", col.names=NA)
#norm gene rpkms
write.table(ALL_rpkm_gene_table_final, paste0(output_folder, "all_samples_gene_RPKMs-normcounts.tab"), quote=FALSE, sep="\t", col.names=NA)
#norm gene tpms log2
write.table(ALL_tpm_gene_table_log2, paste0(output_folder, "all_samples_gene_log2-TPMs-normcounts.tab"), quote=FALSE, sep="\t", col.names=NA)
#norm gene rpkms log2
write.table(ALL_rpkm_gene_table_log2, paste0(output_folder, "all_samples_gene_log2-RPKMs-normcounts.tab"), quote=FALSE, sep="\t", col.names=NA)

