###### Script to save the fine tuned loadings from the spls-DA model.
#Upload libraries
library(mixOmics)
library(hashmap)

#Read arguments
Args = commandArgs(trailingOnly=TRUE);
model_file = (Args[1])
expr_file = (Args[2])
PC_number = (Args[3])
metadata_dir = (Args[4])
all_species = (Args[5])
output_dir = (Args[6])
output_prefix = (Args[7])
output_suffix = (Args[8])
output_suffix_all = (Args[9])

all_species = strsplit(all_species, ",")[[1]]
PC_number = as.numeric(PC_number)

#MAIN
#Load the model
load(model_file)
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]
#I have to find the way to correctly manipulate this
#select.keepX[[1]] = 10

#Get metadata df
all_metadata_df = vector()
for (species in all_species) {
  file_species = species
  if (species == "BmA") {file_species = "Bmo"}
  metadata_df = unique(read.delim(paste0(metadata_dir, "/", file_species, "_samples_info.tab"), header=TRUE)[,1:3])
  colnames(metadata_df) = c("species", "tissues", "metasample")
  metadata_df$species = rep(species, nrow(metadata_df)) #This is necessary only to translate Bmo to BmA
  rownames(metadata_df) = paste0(rep(species, nrow(metadata_df)), "_", as.vector(metadata_df$metasample))
  metadata_df$metasample = NULL
  all_metadata_df = rbind(all_metadata_df, metadata_df)
}
colnames(all_metadata_df) = c("species", "tissues")

#Get expression df
expr_df = read.table(expr_file, header=TRUE, row=1)
expr_df = expr_df[, colnames(expr_df) %in% rownames(all_metadata_df)]

#Generate the needed membership factor
sample_tissue_dict = hashmap(rownames(all_metadata_df), as.vector(all_metadata_df$tissue))
samples_membership = factor(sample_tissue_dict$find(colnames(expr_df)))
splsda_res <- splsda(t(as.matrix(expr_df)), samples_membership, ncomp = ncomp, keepX = select.keepX)

#Select fine-tuned loadings and save to file.
selected_loadings = selectVar(splsda_res, comp = PC_number)$value
write.table(selected_loadings, paste0(output_dir, output_prefix, PC_number, output_suffix), col.names=FALSE, quote=FALSE, sep="\t")
#Select all loadings and save to file.
all_loadings = as.data.frame(splsda_res$loadings.star[[1]][,PC_number])
write.table(all_loadings, paste0(output_dir, output_prefix, PC_number, output_suffix_all), col.names=FALSE, quote=FALSE, sep="\t")

