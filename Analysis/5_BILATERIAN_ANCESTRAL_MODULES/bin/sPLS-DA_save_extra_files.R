#Upload libraries
library(mixOmics)
library(hashmap)

#Read arguments
Args = commandArgs(trailingOnly=TRUE);
model_file = (Args[1])
expr_file = (Args[2])
metadata_dir = (Args[3])
all_species = (Args[4])
output_PC_coords = (Args[5])
output_variance = (Args[6])
output_loadings_number = (Args[7])

all_species = strsplit(all_species, ",")[[1]]

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

PC_df = as.data.frame(splsda_res$variates$X)
explained_variance_df = as.data.frame(splsda_res$prop_expl_var$X)
colnames(PC_df) = sub("comp", "PC", colnames(PC_df))
rownames(explained_variance_df) = sub("comp", "PC", rownames(explained_variance_df))
loadings_number_df = as.data.frame(select.keepX)
rownames(loadings_number_df) = sub("comp", "PC", rownames(loadings_number_df))

#Save files necessary for plotting
write.table(PC_df, output_PC_coords, col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
write.table(explained_variance_df, output_variance, col.names=FALSE, row.names=TRUE, quote=FALSE, sep="\t")
write.table(loadings_number_df, output_loadings_number, col.names=FALSE, row.names=TRUE, quote=FALSE, sep="\t")