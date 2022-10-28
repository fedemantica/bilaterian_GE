#upload libraries
library("ggplot2")
library("mixOmics")
library("hashmap")

#read arguments
Args = commandArgs(trailingOnly=TRUE)
expression_file = (Args[1])
metadata_dir = (Args[2])
all_species = (Args[3])
vertebrata = (Args[4])
insecta = (Args[5])
my_tissue = (Args[6])
output_file = (Args[7])

#upload metadata table
all_species = strsplit(all_species, ",")[[1]]
vertebrata = strsplit(vertebrata, ",")[[1]]
insecta = strsplit(insecta, ",")[[1]]

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
clade_tissue_df = all_metadata_df
clade_tissue_df$clade_tissue = rep("Others", nrow(clade_tissue_df))
clade_tissue_df$clade_tissue[clade_tissue_df$species %in% vertebrata & clade_tissue_df$tissues==my_tissue] = paste0("Vertebrata_", my_tissue)
clade_tissue_df$clade_tissue[clade_tissue_df$species %in% insecta & clade_tissue_df$tissues==my_tissue] = paste0("Insecta_", my_tissue)

sample_tissue_dict = hashmap(rownames(clade_tissue_df), as.vector(clade_tissue_df$clade_tissue))

#upload expression table
expr_df = read.delim(expression_file, header=TRUE, row=1)
#generate a membership factor.
samples_membership = as.factor(sample_tissue_dict$find(colnames(expr_df)))

#list of numbers of features to test for each component
#Setting ncomp to 24, becuase I want results for all categories (Neural_Vertebrata, Neural_Insecta, Neural_Others etc).
list.keepX <- c(1:10, seq(20, 300, 10))
tune.splsda.srbct <- tune.splsda(t(as.matrix(expr_df)), samples_membership, ncomp = 3, validation = 'Mfold', folds = 5, progressBar = TRUE, dist = 'max.dist', measure = "BER", test.keepX = list.keepX, nrepeat=10)

#NB: the output should have a RDATA extension
save(tune.splsda.srbct, file=output_file)
