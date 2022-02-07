##### Script to plot the PCA results and compare between different normalization methods

#Upload libraries
library(ggplot2)
library(impute)
library(cowplot)
library(ggpubr)

#Set arguments
Args = commandArgs(trailingOnly=TRUE);
expr_file = (Args[1])
clade = (Args[2])
all_species = (Args[3])
metadata_dir = (Args[4]) #used to upload the metadata files.
metadata_suffix = (Args[5]) #used to upload the metadata files.
plot_title = (Args[6])
output_file = (Args[7])

#Define functions
PCA_basic_plot = function(expr_df, features_table, species_set, cluster_number, quantity) {
  imputation_out = impute.knn(as.matrix(expr_df), k = 10) #impute missing values.
  imputation_df = as.data.frame(imputation_out$data) #transform to a dataframe.
  #compute PCA. I need the genes on columns in the input table.
  PCA_res = prcomp(as.matrix(t(imputation_df)), center = TRUE, scale = TRUE) #center and scale the PCs.
  ordered_PCA_res = as.data.frame(PCA_res$x[order(rownames(PCA_res$x)),]) 
  #filter the sample only for the samples that made it after the sva correction
  filtered_features_table = features_table[rownames(features_table) %in% colnames(imputation_df),]
  ordered_features_table = filtered_features_table[order(rownames(filtered_features_table)),] #order the table
  all(rownames(ordered_features_table) == rownames(ordered_PCA_res))
  #compute the percentage of variance explained
  eigs <- PCA_res$sdev^2
  first_PC_var = round(eigs[1]/sum(eigs), 3)*100
  second_PC_var = round(eigs[2]/sum(eigs), 3)*100
  third_PC_var = round(eigs[3]/sum(eigs), 3)*100
  fourth_PC_var = round(eigs[4]/sum(eigs), 3)*100
  #filter shapes vector
  my_filtered_shapes_vector=my_shapes_vector[unique(as.vector(ordered_features_table$species))] #filter the shape vector
  #filter ordered species
  filt_ordered_species = ordered_species[ordered_species %in% unique(as.vector(ordered_features_table$species))]
  ordered_features_table$species <- factor(ordered_features_table$species, levels=filt_ordered_species) #change the order of the levels in the factor
  #plot the first two PC
  ordered_PCA_res$tissues=ordered_features_table$tissues
  ordered_PCA_res$species=ordered_features_table$species
  
  p = ggplot() + geom_point(data=as.data.frame(ordered_PCA_res), 
                            aes(x=PC1, y=PC2, 
                                color=tissues, 
                                fill=tissues,
                                shape=species), size=2.5, stroke=0.8) +
    ggtitle(paste0(species_set, " ~ ", cluster_number, ": ", quantity)) + 
    xlab(paste0("PC1: ", first_PC_var, "%")) + 
    ylab(paste0("PC2: ", second_PC_var, "%")) + theme_bw() +
    theme(axis.title = element_text(color="black", size=14), 
          axis.text = element_text(color="black", size=12),
          legend.text = element_text(color="black", size=12), 
          legend.title = element_text(color="black", size=13)) +
    scale_colour_manual(values=my_color_vector) + #this is great because I can have named vectors
    scale_fill_manual(values=my_color_vector) +
    scale_shape_manual(values=my_filtered_shapes_vector) + 
    guides(color=guide_legend(title="Tissues"), shape=guide_legend(title="Species")) 
  #+ coord_fixed()
  
  
  p1 = ggplot() + geom_point(data=as.data.frame(ordered_PCA_res), 
                             aes(x=PC2, y=PC3, 
                                 color=tissues,
                                 fill=tissues,
                                 shape=species), size=2.5, stroke=0.8) +
    ggtitle(paste0(species_set, " ~ ", cluster_number, ": ", quantity)) + 
    xlab(paste0("PC2: ", second_PC_var, "%")) + 
    ylab(paste0("PC3: ", third_PC_var, "%")) + 
    theme_bw() +
    theme(axis.title = element_text(color="black", size=14), 
          axis.text = element_text(color="black", size=12),
          legend.text = element_text(color="black", size=12), 
          legend.title = element_text(color="black", size=13)) +
    scale_colour_manual(values=my_color_vector) + #this is great because I can have named vectors
    scale_fill_manual(values=my_color_vector) +
    scale_shape_manual(values=my_filtered_shapes_vector) + 
    guides(color=guide_legend(title="Tissues"), shape=guide_legend(title="Species"), fill=FALSE)
  
  p2 = ggplot() + geom_point(data=as.data.frame(ordered_PCA_res), 
                             aes(x=PC3, y=PC4, 
                                 color=tissues,
                                 fill=tissues,
                                 shape=species), size=2.5, stroke=0.8) +
    ggtitle(paste0(species_set, " ~ ", cluster_number, ": ", quantity)) + xlab(paste0("PC3: ", third_PC_var, "%")) + ylab(paste0("PC4: ", fourth_PC_var, "%")) + theme_bw() +
    theme(axis.title = element_text(color="black", size=14), 
          axis.text = element_text(color="black", size=12),
          legend.text = element_text(color="black", size=12), 
          legend.title = element_text(color="black", size=13)) +
    scale_colour_manual(values=my_color_vector) + #this is great because I can have named vectors
    scale_fill_manual(values=my_color_vector) +
    scale_shape_manual(values=my_filtered_shapes_vector) + 
    guides(color=guide_legend(title="Tissues"), shape=guide_legend(title="Species"), fill=FALSE)
  
  res = list()
  res[[1]] = p
  res[[2]] = p1
  res[[3]] = p2
  
  my_plots = plot_grid(p, p1, p2, rel_widths=c(1,1,1), align="h", ncol=3)
  return(list("plots" = my_plots, "results"=res))
}

#Set aestetics vectors
my_color_vector = c("Adipose"="firebrick1", "Neural"="mediumblue", "Muscle"="springgreen3", 
                    "Ovary"="darkorchid", "Kidney"="gold", "Testis"="violet", 
                    "DigestiveTract"="chocolate1", "Epithelial"="steelblue1")
my_shapes_vector = c("Hs2"=0, "Mm2"=1, "Bt2"=2, "Mdo"=3, "Gga"=4, "Xtr"=5, "Dre"=6, "Cmi"=7, "Bla"=9, "Sp2"=8,
                     "Dme"=15, "Eba"=18, "Aae"= 23, "BmA"=25, "Tca"=17, "Bge"=10, "Ame"=16, "Cdi"=13, "Sma"=14, "Obi"=11)
my_alpha_vector=c("Hs2"=1, "Mm2"=1, "Bt2"=1, "Mdo"=1, "Gga"=1, "Xtr"=1, "Dr2"=1, "Cmi"=1, "Bla"=0.5, "Spu"=0.5, 
                  "Sma"=0.5, "Ame"=0.5, "Dme"=0.5, "Cdi"=0.5, "Obi"=0.5, "Bge"=0.5, "BmA"=0.5, "Tca"=0.5)
my_sets_color = c("Vertebrata"= "mediumorchid", "Deuterostoma"="olivedrab3", "Insecta"="tan3", "Protostoma"="yellow", "Bilateria"="firebrick1")
my_species_colors=c("Hs2"="#00BFFF", "Mm2"="#159DFF", "Bt2"="#1873CD", "Mdo"="#174B8B", "Gga"="#3F3F8B", 
                    "Xtr"="#64398B", "Dre"="#82359D", "Cmi"="#9932CC", "Bla"="#2E8B57", "Sp2"="#96B12B", 
                    "Obi"="#FFFF00", "Sma"="#FFCC33", "Cdi"="goldenrod", "Bge"="burlywood2", "Ame"="#FF9966", 
                    "Tca"="#FF6600", "BmA"="#CC6600", "Aae"="#993300", "Eba"="firebrick3", "Dme"="red")

#Read the input file
expr_df = read.delim(expr_file, sep="\t", header=TRUE, row=1)
expr_df = expr_df[complete.cases(expr_df),]

#Define variables
all_species = strsplit(all_species, ",")[[1]] #Get vector from string of comma separated elements
ordered_species = all_species

#Generate metadata
all_metadata_df = vector()
for (species in all_species) {
  file_species = species
  if (species == "BmA") {file_species= "Bmo"}
  metadata_df = unique(read.delim(paste0(metadata_dir, file_species, metadata_suffix), header=TRUE)[,1:3])
  colnames(metadata_df) = c("species", "tissues", "metasample")
  metadata_df$species = rep(species, nrow(metadata_df)) #This is necessary only to translate Bmo to BmA
  rownames(metadata_df) = paste0(rep(species, nrow(metadata_df)), "_", as.vector(metadata_df$metasample))
  metadata_df$metasample = NULL
  all_metadata_df = rbind(all_metadata_df, metadata_df)
}
colnames(all_metadata_df) = c("species", "tissues")

#Run PCA
orthogroups_number=nrow(expr_df)
res_PCA = PCA_basic_plot(expr_df, all_metadata_df, clade, orthogroups_number, plot_title) #This plot title is a useless argument. Just don't feel like cleaning up now.

#Plot
PCA_plot = annotate_figure(ggarrange(res_PCA$results[[1]] + ggtitle(""), 
                                      res_PCA$results[[2]] + ggtitle(""), 
                                      common.legend=TRUE, legend="right", 
                                      ncol=2), 
                            fig.lab=paste0(plot_title, " ~ ", orthogroups_number, " OGs"), fig.lab.pos = "top.left", fig.lab.size=18)

#Save plot to output file
pdf(output_file, width = 12, height = 7)
PCA_plot
dev.off()
#ggsave(PCA_plot, output_file, "pdf")