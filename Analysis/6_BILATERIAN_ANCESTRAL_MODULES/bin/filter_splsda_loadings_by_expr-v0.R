### This is a script to select those loadings which present the highest expression in the tissue of interest
######## Upload libraries
library("reshape2")
library("tidyverse")

######## Set Arguments
Args <- commandArgs(trailingOnly=TRUE);
expr_file = (Args[1]);
tissues = (Args[2]);
vertebrata = (Args[3]);
insecta = (Args[4]);
outgroups = (Args[5]);
loadings_dir = (Args[6]);
loadings_suffix = (Args[7]);
output_dir = (Args[8]);
output_suffix_all = (Args[9]);
output_suffix_Hs2_GO_input = (Args[10]);
output_suffix_orthogroups_GO_input = (Args[11]);

#####This is an inside variable
tissue_PC_vector = c("Neural"="Negative_PC1",
                     "Testis"="Positive_PC2",
                     "Ovary"="Negative_PC3",
                     "Muscle"="Positive_PC4",
                     "DigestiveTract"="Negative_PC5",
                     "Epithelial"="Positive_PC6",
                     "Kidney"="Negative_PC7",
                     "Adipose"="Negative_PC8")

#tissue_PC_vector = c("Neural"="Negative_PC1",
#                     "Testis"="Positive_PC2",
#                     "Ovary"="Negative_PC4",
#                     "Muscle"="Positive_PC3",
#                     "DigestiveTract"="Negative_PC6",
#                     "Epithelial"="Positive_PC5",
#                     "Kidney"="Positive_PC6",
#                     "Adipose"="Negative_PC7")

#Generate variable vectors from strings
tissues = strsplit(tissues, ",")[[1]]
vertebrata = strsplit(vertebrata, ",")[[1]]
insecta = strsplit(insecta, ",")[[1]]
outgroups = strsplit(outgroups, ",")[[1]]

######### Read input table
expr_df = read.table(expr_file, header=TRUE, row=1)

#Transform input table to long format
expr_df_long = melt(as.matrix(expr_df))
colnames(expr_df_long) = c("OG_ID", "Species_Tissue", "Expr")
expr_df_long$Species = sub("_.*", "", expr_df_long$Species_Tissue)
expr_df_long$Tissue = sub(".*_", "", expr_df_long$Species_Tissue)

#Add clade so that the dataframe could be grouped by that.
clade_expr_df = expr_df_long %>% 
  mutate(Clade = case_when(Species %in% vertebrata ~ "Vertebrata",
                           Species %in% insecta ~ "Insecta",
                           Species %in% outgroups ~ "Outgroups")
  )

######## Compute zscore
#Compute the zscore separately for Insects, Vertebrate and Outgroups across tissues
OG_zscore_df = clade_expr_df %>%
  group_by(OG_ID, Tissue, Clade) %>%
  summarize(zscore=median(Expr, na.rm=TRUE)) %>%
  ungroup %>%
  #group_by(OG_ID, Clade) %>% 
  #mutate(zscore = (median_expr - mean(median_expr, na.rm=TRUE))/sd(median_expr, na.rm=TRUE)) %>%
  unite("OG_ID_Tissue", OG_ID, Tissue, remove=FALSE, sep="-")

#For each tissue, select all those genes where the tissue of interest has the highest Zscore in Vertebrate and Insects.
tissue_maximum_all_clades_raw_df = OG_zscore_df %>%
  group_by(OG_ID, Clade) %>%
  filter(zscore == max(zscore)) %>%
  arrange(OG_ID, Clade) %>%
  ungroup

#Subset to only Vertebrates and Insects highest value in the tissue
tissue_maximum_two_clades_raw_df = subset(tissue_maximum_all_clades_raw_df, Clade %in% c("Vertebrata", "Insecta"))
tissue_maximum_all_clades = names(table(tissue_maximum_two_clades_raw_df$OG_ID_Tissue)[table(tissue_maximum_two_clades_raw_df$OG_ID_Tissue)==2])
tissue_maximum_all_clades_df = subset(tissue_maximum_all_clades_raw_df, OG_ID_Tissue %in% tissue_maximum_all_clades)


###################################
####### Save loadings #############
###################################

for (query_tissue in tissues) {
  #### Upload spls-DA loadings
  file_prefix = unname(tissue_PC_vector[query_tissue])
  my_loadings_df = read.delim(paste0(loadings_dir, "/", file_prefix, loadings_suffix), 
                              header=FALSE, col.names=c("OG_ID", "Hs2_ID", "Dme_ID", "Mm2_ID"))
  loadings_OGs = as.vector(my_loadings_df$OG_ID)
  
  #### Get gene orthogroups with the highest zscore in each tissue
  tissue_zscore_selected_df = subset(tissue_maximum_all_clades_df, Tissue==query_tissue)
  tissue_zscore_selected_OGs = unique(as.vector(tissue_zscore_selected_df$OG_ID))
  
  #Count OGs retrieved by both approaches
  common_OGs_zscore = loadings_OGs[loadings_OGs %in% tissue_zscore_selected_OGs]
  #Get the corresponding human gene
  common_human_genes_zscore = as.vector(subset(my_loadings_df, OG_ID %in% common_OGs_zscore)$Hs2_ID)
  common_human_genes_zscore = common_human_genes_zscore[!(is.na(common_human_genes_zscore))]
  #common_human_genes_zscore = vector()
  #for (OG_ID in common_OGs_zscore) {
    #human_gene = geneID_names_dict$keys()[geneID_names_dict$values()==OGID_humanID_dict$find(OG_ID)]
    #common_human_genes_zscore = c(common_human_genes_zscore, human_gene)
  #}
  #Save orthogroup IDs and gene IDs to file
  final_loadings_df = subset(my_loadings_df, OG_ID %in% common_OGs_zscore)
  write.table(final_loadings_df, paste0(output_dir, "/", query_tissue, output_suffix_all),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  #Save human geneID to file
  write.table(common_human_genes_zscore, paste0(output_dir, "/", query_tissue, output_suffix_Hs2_GO_input),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  #Save orthogroups IDs to file
  write.table(common_OGs_zscore, paste0(output_dir, "/", query_tissue, output_suffix_orthogroups_GO_input),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
