#upload libraries
library("ggplot2")
library("mixOmics")
library("hashmap")

#read arguments
Args = commandArgs(trailingOnly=TRUE)
expression_file = (Args[1])
metadata_dir = (Args[2])
clade_species = (Args[3])
output_file = (Args[4])

#upload metadata table
all_species = strsplit(clade_species, ",")[[1]]

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
#generate sample-tissue dictionary
sample_tissue_dict = hashmap(rownames(all_metadata_df), as.vector(all_metadata_df$tissues))

#upload expression table
expr_df = read.delim(expression_file, header=TRUE, row=1)
#generate a membership factor. Given the results of the PCA, I group Testis and Ovary together.
samples_membership = sample_tissue_dict$find(colnames(expr_df))

#list of numbers of features to test for each component
list.keepX <- c(1:10, seq(20, 300, 10))
tune.splsda.srbct <- tune.splsda(t(as.matrix(expr_df)), samples_membership, ncomp = 10, validation = 'Mfold', folds = 4, progressBar = TRUE, dist = 'max.dist', measure = "BER", test.keepX = list.keepX, nrepeat=10)

#NB: the output should have a RDATA extension
save(tune.splsda.srbct, file=output_file)
