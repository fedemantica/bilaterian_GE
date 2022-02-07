##Script to the trees from which we start building the branches lenghts
#Upload libraries
library(ggplot2)
library(ggtree)

#Read arguments
Args = commandArgs(trailingOnly=TRUE);
expr_file = (Args[1])
tissues = (Args[2])
output_file = (Args[3])

#Modify arguments
ordered_tissues = strsplit(tissues, ",")[[1]]

#Define functions
number_ticks <- function(n) {function(limits) pretty(limits, n)}

get_median_tissue_tree = function(sample_expr_table, species_set, my_tissue, n_bootstraps) {
  #select list of Neural samples for all the represented species
  set_species = unique(gsub("_.*", "", colnames(sample_expr_table)))
  #select the tissue columns belonging to the set species
  tissue_samples = colnames(sample_expr_table)[grep(my_tissue, colnames(sample_expr_table))]
  #subset the table: only tissue samples
  tissue_expr_table = sample_expr_table[,tissue_samples]
  #compute distance matrix from the transposed  dataframe
  tissue_distance_matrix = dist(as.matrix(t(tissue_expr_table)), method="euclidean")
  #get NJ tree
  tissue_tree = nj(tissue_distance_matrix)
  
  #perform 1000 bootstraps
  bootstrap = boot.phylo(tissue_tree, as.matrix(t(tissue_expr_table)), function(x) nj(dist(x)), n_bootstraps)
  corrected_bootstrap = bootstrap/n_bootstraps
  
  my_tree = ggtree(tissue_tree) +
    ggtitle(paste0(my_tissue, " ~ ", species_set)) +
    geom_tiplab() + geom_treescale() + geom_nodelab(label=corrected_bootstrap, geom="text", hjust=-.2, color="coral3", face="bold")
  return(my_tree)
}

###################################
########## MAIN ###################
###################################

#Read inputs

#Apply function
expr_table = read.delim(expr_file, header=TRUE, row=1)
for (tissue in ordered_tissues) {
  tree_plot = get_median_tissue_tree(expr_table, species_set, tissue, 1000)
  print(tree_plot)
}

#Save to output
pdf(output_file)
print(tree_plot)
dev.off()