######## Script to plot the tissue expression correlations across species of the specified clade

#Upload libraries
library(ggplot2)
library(reshape2)
library(gridExtra)

#Read arguments
Args = commandArgs(trailingOnly=TRUE);
expr_file = (Args[1])
query_species = (Args[2])
tissues = (Args[3])
species = (Args[4])
output_file = (Args[5])

#Modify arguments
ordered_tissues = strsplit(tissues, ",")[[1]]
ordered_species = strsplit(species, ",")[[1]]

#Define functions
number_ticks <- function(n) {function(limits) pretty(limits, n)}

get_cor_table = function(expr_table, query_species, tissue)  {
  #select the Neural columns belonging to the set species
  tissue_samples = colnames(expr_table)[grep(tissue, colnames(expr_table))]
  #subset the table: only tissue samples
  tissue_expr_table = expr_table[,tissue_samples]
  colnames(tissue_expr_table) = sub("_.*", "", colnames(tissue_expr_table))
  clade_species = colnames(tissue_expr_table)
  #get target species depending on the clade
  all_target_species = clade_species
  final_matrix =  vector()
  if (!(query_species %in% colnames(tissue_expr_table))) {
    final_table = c(query_species, NA, tissue, NA)
    return(final_table)} else {
      for (target_species in all_target_species) {
        subsetted_table = tissue_expr_table[,c(query_species, target_species)]
        subsetted_table = subsetted_table[complete.cases(subsetted_table),] #remove all rows containing one or two NAs
        query_vector  = as.vector(subsetted_table[,query_species]) #isolate query vector
        target_vector = as.vector(subsetted_table[,target_species]) #isolate target vector
        rho = cor.test(query_vector, target_vector, method="spearman")$estimate #isolate rho
        final_matrix = rbind(final_matrix, c(query_species, target_species, tissue, rho)) #append to final matrix
      }
      final_table = as.data.frame(final_matrix)
      colnames(final_table) = c("species_query", "species_target", "tissue", "cor")
      final_table$cor = as.numeric(as.vector(final_table$cor))
      return(final_table)
    }
}

#Set aestetics vectors
tissue_colors_vector = c("Adipose"="firebrick1", "Neural"="mediumblue", "Muscle"="springgreen3", 
                    "Ovary"="darkorchid", "Kidney"="gold", "Testis"="violet", 
                    "DigestiveTract"="chocolate1", "Epithelial"="steelblue1")
species_shapes_vector = c("Hs2"=0, "Mm2"=1, "Bt2"=2, "Mdo"=3, "Gga"=4, "Xtr"=5, "Dre"=6, "Cmi"=7, "Bla"=9, "Sp2"=8,
                     "Dme"=15, "Eba"=18, "Aae"= 23, "BmA"=25, "Tca"=17, "Bge"=10, "Ame"=16, "Cdi"=13, "Sma"=14, "Obi"=11)
species_alphas_vector = c("Hs2"=1, "Mm2"=1, "Bt2"=1, "Mdo"=1, "Gga"=1, "Xtr"=1, "Dr2"=1, "Cmi"=1, "Bla"=0.5, "Spu"=0.5, 
                  "Sma"=0.5, "Ame"=0.5, "Dme"=0.5, "Cdi"=0.5, "Obi"=0.5, "Bge"=0.5, "BmA"=0.5, "Tca"=0.5)
clade_colors_vector = c("Vertebrata"= "mediumorchid", "Deuterostoma"="olivedrab3", "Insecta"="tan3", "Protostoma"="yellow", "Bilateria"="firebrick1")
species_colors_vector = c("Hs2"="#00BFFF", "Mm2"="#159DFF", "Bt2"="#1873CD", "Mdo"="#174B8B", "Gga"="#3F3F8B", 
                    "Xtr"="#64398B", "Dre"="#82359D", "Cmi"="#9932CC", "Bla"="#2E8B57", "Sp2"="#96B12B", 
                    "Obi"="#FFFF00", "Sma"="#FFCC33", "Cdi"="goldenrod", "Bge"="burlywood2", "Ame"="#FF9966", 
                    "Tca"="#FF6600", "BmA"="#CC6600", "Aae"="#993300", "Eba"="firebrick3", "Dme"="red")

################################
###### Main ####################
################################
#Upload expression table
expr_full_table = read.table(expr_file, header=TRUE, row=1)
expr_table = expr_full_table[,grep(paste0(ordered_species, collapse="|"), colnames(expr_full_table), value=TRUE)] #subset in case the input is for all Bilateria

query_species_table = vector()
for (tissue in ordered_tissues) {
  tissue_cor_table = get_cor_table(expr_table, query_species, tissue)
  query_species_table = rbind(query_species_table, tissue_cor_table)
}
query_species_table$species_target = factor(query_species_table$species_target, levels = ordered_species)

###### Query species plot #########
query_species_cor_plot = ggplot() +
  geom_point(data=query_species_table, aes(x=species_target, y=cor, color=tissue)) +
  geom_line(data=query_species_table, aes(x=species_target, y=cor, color=tissue, group=tissue)) +
  theme_bw() +
  theme(axis.text.x =element_text(color="black", size=12, angle=90, hjust=1, vjust=0.5),
        axis.text.y =element_text(color="black", size=12)) +
  scale_color_manual(values=tissue_colors_vector) +
  ggtitle(paste0(query_species, " as query")) +
  scale_y_continuous(breaks=number_ticks(10), limits=c(0,1)) +
  guides(color=FALSE)

##### All species plot ############

all_species_table = vector()
for (tissue in ordered_tissues) {
  for (new_query_species in ordered_species) {
    tissue_cor_table = get_cor_table(expr_table, new_query_species, tissue)
    all_species_table = rbind(all_species_table, tissue_cor_table)
  }
}

#Subset to just one value per species combination (currently: 2)
all_species_table$species_comb = paste0(all_species_table$species_query, "_", all_species_table$species_target)
all_species_comb = paste0(as.data.frame(t(combn(ordered_species, 2)))[,1], "_", as.data.frame(t(combn(ordered_species, 2)))[,2]) #all combinations between species
all_species_table = subset(all_species_table, species_comb %in% all_species_comb)
all_species_table$species_comb = factor(all_species_table$species_comb, levels = all_species_comb) #Order factor
all_species_table$cor = as.numeric(as.vector(all_species_table$cor)) #Transform correlations to a numeric vector

### Boxplot: x=tissue, y=correlations (each dot corresponds to a species pair)
plot_list = list()
all_species_cor_boxplot = ggplot(data=all_species_table, aes(x=tissue, y=cor, fill=tissue)) +
  geom_jitter(color="black", size=0.3) +
  geom_boxplot(alpha=0.7, color="black") +
  theme_bw() +
  scale_fill_manual(values=tissue_colors_vector) +
  theme(axis.text.x=element_text(color="black", angle=30, hjust=1, vjust=1),
        axis.text.y=element_text(color="black")) +
  guides(fill=FALSE) +
  ylim(c(0.40,0.80)) +
  ggtitle("All species pairs ~ Spearman's correlations by tissue")

#add to plot_list
plot_list[[1]] = all_species_cor_boxplot

### Barplot: x=Species pair, y=correlation (one panel per tissue)
i = 2
for (tissue in ordered_tissues) {
  tissue_table = subset(all_species_table, tissue == tissue)
  all_species_cor_barplot = ggplot(data=tissue_table, aes(x=species_comb, y=cor, fill=tissue), color="white") +
    geom_bar(alpha=0.8, stat="identity") + 
    scale_fill_manual(values=tissue_colors_vector) +
    theme_bw() +
    theme(axis.text.x=element_text(color="black", angle=45, hjust=1, vjust=1),
          axis.text.y=element_text(color="black")) +
    ggtitle(paste0(tissue, " ~ Spearman's correlations")) +
    guides(fill=FALSE) +
    ylim(c(0,1))
  
  plot_list[[i]] = all_species_cor_barplot
  i = i+1
}

#### Save to output file #########
pdf(output_file)
print(query_species_cor_plot)
do.call("grid.arrange", c(plot_list, ncol=3))
dev.off()