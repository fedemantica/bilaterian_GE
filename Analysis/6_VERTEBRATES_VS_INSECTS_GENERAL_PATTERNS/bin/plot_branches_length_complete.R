######## Script to plot the branches lengths analysis results

#Upload libraries
library(ggplot2)
library(ape)

#Read arguments
Args = commandArgs(trailingOnly=TRUE);
expr_file = (Args[1])
clade = (Args[2])
tissues = (Args[3])
species = (Args[4])
output_lengths_file = (Args[5])
output_plot_file = (Args[6])

#Modify arguments
ordered_tissues = strsplit(tissues, ",")[[1]]
ordered_species = strsplit(species, ",")[[1]]
collapsed_species = gsub(",","|",species)

#Define functions
number_ticks <- function(n) {function(limits) pretty(limits, n)}

extract_branch_lengths = function(tissue_expr_table, collapsed_species, tissue, n_bootstraps) {
  #subset to only species in the set
  tissue_expr_table = tissue_expr_table[,grep(collapsed_species, colnames(tissue_expr_table), value=TRUE)]
  #select the tissue columns belonging to the set species
  tissue_samples = colnames(tissue_expr_table)[grep(tissue, colnames(tissue_expr_table))]
  #subset the table: only neural samples
  tissue_expr_table = tissue_expr_table[,tissue_samples]
  if (exists("tissue_expr_table")) {
    #compute distance matrix from the transposed  dataframe
    tissue_distance_matrix = cor(tissue_expr_table, method="spearman", use="pairwise.complete.obs")
    tissue_distance_matrix = apply(tissue_distance_matrix, c(1,2), function(x) 1-x)
    #get NJ tree
    tissue_tree = nj(tissue_distance_matrix)
    #perform 1000 bootstraps
    bootstrap = boot.phylo(tissue_tree, as.matrix(t(tissue_expr_table)), function(x) nj(dist(x)), n_bootstraps, trees=TRUE)
    branches_lengths = unlist(lapply(bootstrap$trees, function(x) sum(x$edge.length)))
    return(branches_lengths)
  }
  else {return("not_enough_samples")}
}

#Set aestetics vectors
tissue_colors_vector = c("Adipose"="firebrick1", "Neural"="mediumblue", "Muscle"="springgreen3", 
                         "Ovary"="darkorchid", "Kidney"="gold", "Testis"="violet", 
                         "DigestiveTract"="chocolate1", "Epithelial"="steelblue1")

################################
###### Main ####################
################################
tissue_expr_table = read.delim(expr_file, header=TRUE, row=1) #Read expression table in.
tissue_expr_table = tissue_expr_table[complete.cases(tissue_expr_table),]

##### Save branches lengths #######

#Open output file
pdf(output_plot_file)

for (subclade in c("Vertebrata", "Insecta")) {
  clade_values_df = vector()
  if (subclade=="Vertebrata") {collapsed_species = c("Hs2|Mm2|Bt2|Mdo|Gga|Xtr|Dre|Cmi")}
  if (subclade=="Insecta") {collapsed_species = c("Dme|Eba|Aae|Bmo|Tca|Ame|Bge|Cdi")}
  
  for (tissue in ordered_tissues) {
    all_tree_lengths = extract_branch_lengths(tissue_expr_table, collapsed_species, tissue, 1000)
    all_tree_lengths_df = data.frame(set=rep(clade, length(all_tree_lengths)), 
                                     tissue=rep(tissue, length(all_tree_lengths)), 
                                     values=all_tree_lengths)
    clade_values_df = rbind(clade_values_df, all_tree_lengths_df)
  }
  #Save table to output file
  write.table(clade_values_df, output_lengths_file, sep="\t", quote=FALSE, row.names = FALSE)
  
  ##### Plot branches lengths #######
  #Reformat input df (just to be sure) 
  clade_values_df$tissue = factor(clade_values_df$tissue, levels=ordered_tissues)
  clade_values_df$values = as.numeric(as.vector(clade_values_df$values))
  
  ##### Plot
  branches_length_plot = ggplot() +
    geom_boxplot(data=clade_values_df, aes(x=tissue, y=values, fill=tissue, color=tissue), alpha=0.7, outlier.shape = NA) +
    theme_bw() +
    scale_fill_manual(values=tissue_colors_vector) +
    scale_color_manual(values=tissue_colors_vector) +
    ggtitle(clade) +
    theme(axis.text.x = element_text(color="black", size=12, angle=30, hjust = 1),
          axis.text.y = element_text(color="black", size=12),
          axis.title = element_text(color="black", size=14),
          plot.title = element_text(color="black", size=16, hjust=0.5)) +
    ylab("Tot branch length") +
    xlab("Tissue")  +
    guides(fill=FALSE,  color=FALSE) +
    scale_y_continuous(breaks=number_ticks(8))
  
  print(branches_length_plot)
}

dev.off()
