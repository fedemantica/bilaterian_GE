####### Rscript to get normalized counts

###### Upload libraries
library("DESeq2")
library("ggplot2")

#### Read arguments
Args <- commandArgs(trailingOnly=TRUE);
counts_file = Args[1]
metadata_dir = Args[2]
tissue = Args[3]
species_string = Args[4]
output_dir = Args[5]

#### Generate variables
all_species = strsplit(species_string, ",")[[1]]

#### Read counts file
input_df = read.delim(counts_file, header=TRUE, row=1)
input_matrix = as.matrix(input_df)
input_matrix_rounded = apply(input_matrix, c(1,2), round) #Round counts that come from an average.
input_matrix_complete = input_matrix_rounded[complete.cases(input_matrix_rounded),] #select only the gene orthogroups where there are no NAs.

##### Generate metadata
all_metadata_df = vector()
for (species in all_species) {
  file_species = species
  if (species == "BmA") {file_species= "Bmo"} #This is just because I was using two different species identifiers
  metadata_df = unique(read.delim(paste0(metadata_dir, file_species, "_samples_info.tab"), header=TRUE)[,1:3])
  colnames(metadata_df) = c("species", "tissues", "metasample")
  metadata_df$species = rep(species, nrow(metadata_df)) #This is necessary only to translate Bmo to BmA
  rownames(metadata_df) = paste0(rep(species, nrow(metadata_df)), "_", as.vector(metadata_df$metasample))
  metadata_df$metasample = NULL
  all_metadata_df = rbind(all_metadata_df, metadata_df)
}
colnames(all_metadata_df) = c("species", "tissues")
all_metadata_df = all_metadata_df[colnames(input_df),]

#### Add column relative to the tissue of interest (my_tissue)
all_metadata_df$my_tissue = rep(tissue, nrow(all_metadata_df))
all_metadata_df$my_tissue[all_metadata_df$tissues != tissue] = paste0("Not_", tissue)


###### run DESeq2
dds <- DESeqDataSetFromMatrix(countData = input_matrix_complete, colData = all_metadata_df, design = ~my_tissue)
dds <- estimateSizeFactors(dds)
#set the other tissues as reference levels
dds$my_tissue = relevel(dds$my_tissue, paste0("Not_", tissue))
#Get the log2FC between condition of interest
dds <- DESeq(dds)
#get Results
my_res = results(dds, alpha = 0.05)

### Build contrast
my_contrast = paste0("my_tissue_", tissue, "_vs_Not_", tissue)
my_res_shr <- lfcShrink(dds, coef=my_contrast, res=my_res)
### Convert results to dataframes
my_res_df = as.data.frame(my_res)
my_res_shr_df = as.data.frame(my_res_shr)

#oreder by pvalue
my_res_shr_df = my_res_shr_df[order(my_res_shr_df$padj),]
#get significantly upregulated in tissue
upregulated_df = subset(my_res_shr_df, padj <= 0.05 & log2FoldChange > 0)
downregulated_df = subset(my_res_shr_df, padj <= 0.05 & log2FoldChange < 0)


### Save dataframe to output files
write.table(my_res_shr_df, paste0(output_dir, "/", tissue, "-All_shrunk_LFC.tab"), sep="\t", quote=F, col.names=NA)
write.table(upregulated_df, paste0(output_dir, "/", tissue, "-UPregulated_shrunk_LFC-All.tab"), sep="\t", quote=F, col.names=NA)
write.table(downregulated_df, paste0(output_dir, "/", tissue, "-DOWNregulated_shrunk_LFC-All.tab"), sep="\t", quote=F, col.names=NA)

## #save the MA plots
pdf(paste0(output_dir, "/", tissue, "-MA_plots.pdf"))
plotMA(my_res, ylim=c(-6,6))
plotMA(my_res_shr, ylim=c(-6,6))
#Add pvalues plots
ggplot() + 
  geom_histogram(data=my_res_shr_df, aes(x=pvalue), binwidth = 0.05, color="white", fill="steelblue") +
  theme_bw()
#close output devide
dev.off()
