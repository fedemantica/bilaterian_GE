###### Script to save the fine tuned loadings from the spls-DA model.
#Upload libraries
library(mixOmics)
library(hashmap)

#Read arguments
Args = commandArgs(trailingOnly=TRUE);
model_file = (Args[1])
expr_file = (Args[2])
tissue = (Args[3])
PC_number = (Args[4])
metadata_dir = (Args[5])
all_species = (Args[6])
vertebrata = (Args[7])
insecta = (Args[8])
output_dir = (Args[9])
output_prefix = (Args[10])
output_suffix = (Args[11])
output_suffix_all = (Args[12])

all_species = strsplit(all_species, ",")[[1]]
vertebrata = strsplit(vertebrata, ",")[[1]]
insecta = strsplit(insecta, ",")[[1]]
PC_number = as.numeric(PC_number)

#MAIN
#Load the model
load(model_file)
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]
#I have to find the way to correctly manipulate this
#Try to increase the number of genes
if (tissue=="Neural") {
  select.keepX[[2]] = 50
}

#for (i in seq(1,ncomp)) {
#  if (select.keepX[[i]] <= 20) {select.keepX[[i]] = 40}
#}

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

#This part is for when Vertebrates and Insects are separated
all_metadata_df$clade_tissue = rep("Others", nrow(all_metadata_df))
all_metadata_df$clade_tissue[all_metadata_df$species %in% vertebrata & all_metadata_df$tissues==tissue] = paste0("Vertebrata_", tissue)
all_metadata_df$clade_tissue[all_metadata_df$species %in% insecta & all_metadata_df$tissues==tissue] = paste0("Insecta_", tissue)

#Generate the needed membership factor
sample_tissue_dict = hashmap(rownames(all_metadata_df), as.vector(all_metadata_df$clade_tissue))
samples_membership = factor(sample_tissue_dict$find(colnames(expr_df)))
splsda_res <- splsda(t(as.matrix(expr_df)), samples_membership, ncomp = ncomp, keepX = select.keepX)

#Select fine-tuned loadings and save to file.
selected_loadings = selectVar(splsda_res, comp = PC_number)$value
write.table(selected_loadings, paste0(output_dir, output_prefix, PC_number, output_suffix), col.names=FALSE, quote=FALSE, sep="\t")
#Select all loadings and save to file.
all_loadings = as.data.frame(splsda_res$loadings.star[[1]][,PC_number])
write.table(all_loadings, paste0(output_dir, output_prefix, PC_number, output_suffix_all), col.names=FALSE, quote=FALSE, sep="\t")
