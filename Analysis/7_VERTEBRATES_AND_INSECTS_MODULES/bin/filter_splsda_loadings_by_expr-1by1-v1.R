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
#And this has to be re-defined for every run
clade_tissue_combs = c("Vertebrata_Neural", "Insecta_Neural",
                       "Vertebrata_Testis", "Insecta_Testis",
                       "Vertebrata_Ovary", "Insecta_Ovary",
                       "Vertebrata_Muscle", "Insecta_Muscle",
                       "Vertebrata_Kidney", "Insecta_Kidney",
                       "Vertebrata_Epithelial", "Insecta_Epithelial",
                       "Vertebrata_DigestiveTract", "Insecta_DigestiveTract",
                       "Vertebrata_Adipose", "Insecta_Adipose")

tissue_PC_list = list("Vertebrata_Neural" = "Positive_PC2",
                     "Insecta_Neural" = "Negative_PC2",
                     "Vertebrata_Testis" = c("Positive_PC2", "Positive_PC3"),
                     "Insecta_Testis" = "Negative_PC2",
                     "Vertebrata_Ovary" = "Positive_PC2",
                     "Insecta_Ovary" = "Negative_PC1",
                     "Vertebrata_Muscle" = "Negative_PC2",
                     "Insecta_Muscle" = c("Positive_PC2", "Positive_PC3"),
                     "Vertebrata_Kidney" = c("Positive_PC1", "Positive_PC3"),
                     "Insecta_Kidney" = "Negative_PC2",
                     "Vertebrata_Epithelial" = "Positive_PC1",
                     "Insecta_Epithelial" = "Negative_PC2",
                     "Vertebrata_DigestiveTract" = "Positive_PC1",
                     "Insecta_DigestiveTract" = "Negative_PC2",
                     "Vertebrata_Adipose" = "Positive_PC1",
                     "Insecta_Adipose" = "Negative_PC2")

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
  unite("OG_ID_Tissue", OG_ID, Tissue, remove=FALSE, sep="-")

#For each tissue and each clade, select all those genes where the tissue of interest has the highest Zscore.
tissue_maximum_all_clades_raw_df = OG_zscore_df %>%
  group_by(OG_ID, Clade) %>%
  filter(zscore == max(zscore)) %>%
  arrange(OG_ID, Clade) %>%
  ungroup

tissue_maximum_two_clades_raw_df = subset(tissue_maximum_all_clades_raw_df, Clade %in% c("Vertebrata", "Insecta"))


###################################
####### Save loadings #############
###################################

######## Cycle on each combination of clade and tissue
for (comb in clade_tissue_combs) {
  my_clade = sub("_.*", "", comb)
  my_tissue = sub(".*_", "", comb)
  input_prefix_vector = unname(tissue_PC_list[comb][[1]])
  loadings_OGs = vector()
  all_loadings_df = vector()

  #Require that the other Clade does not have the highest zscore in the same tissue
  if (my_clade=="Vertebrata") {other_clade="Insecta"}
  if (my_clade=="Insecta") {other_clade="Vertebrata"}
  other_clade_highest_in_tissue_OGs = as.vector(subset(tissue_maximum_all_clades_raw_df, Tissue==my_tissue & Clade==other_clade)$OG_ID)
  #### Get gene orthogroups with the highest zscore in each tissue
  tissue_zscore_selected_df = subset(tissue_maximum_all_clades_raw_df, Tissue==my_tissue & Clade==my_clade)
  tissue_zscore_selected_OGs = unique(as.vector(tissue_zscore_selected_df$OG_ID))
  tissue_zscore_selected_OGs = tissue_zscore_selected_OGs[!(tissue_zscore_selected_OGs %in% other_clade_highest_in_tissue_OGs)]
  #Require that the zscore in the tissue is higher than the zscore in the tissue of the other clade
  tissue_highest_by_clade_df = OG_zscore_df %>%
    group_by(OG_ID, Clade) %>%
    filter(zscore == max(zscore)) %>%
    ungroup %>%
    filter(Clade %in% c("Vertebrata", "Insecta") & Tissue==my_tissue) %>%
    group_by(OG_ID) %>%
    filter(zscore == max(zscore)) %>%
    ungroup %>%
    filter(Clade==my_clade)
  tissue_highest_by_clade_OGs = tissue_highest_by_clade_df$OG_ID
  
  tissue_zscore_selected_OGs = tissue_zscore_selected_OGs[tissue_zscore_selected_OGs %in% tissue_highest_by_clade_OGs]
  
  #for each element in the file prefix vector
  for (input_prefix in input_prefix_vector) {
    my_loadings_df = read.delim(paste0(loadings_dir, "/", my_tissue, "/", input_prefix, loadings_suffix), 
                                header=FALSE, col.names=c("OG_ID", "Hs2_ID", "Dme_ID", "Mm2_ID"))
    loadings_OGs = c(loadings_OGs, my_loadings_df$OG_ID)
    all_loadings_df = rbind(all_loadings_df, my_loadings_df)
  }
  #Select OGs retrieved by both approaches
  common_OGs_zscore = loadings_OGs[loadings_OGs %in% tissue_zscore_selected_OGs]
  #Get the corresponding human gene
  common_human_genes_zscore = as.vector(subset(all_loadings_df, OG_ID %in% common_OGs_zscore)$Hs2_ID)
  common_human_genes_zscore = common_human_genes_zscore[!(is.na(common_human_genes_zscore))]
  #Get final loadings
  final_loadings_df = subset(my_loadings_df, OG_ID %in% common_OGs_zscore)
  #Save orthogroup IDs and gene IDs to file
  write.table(final_loadings_df, paste0(output_dir, "/", my_tissue, "/splsda_zscore_comb/", comb, output_suffix_all),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  #Save human geneID to file
  write.table(common_human_genes_zscore, paste0(output_dir, "/", my_tissue, "/splsda_zscore_comb/", comb, output_suffix_Hs2_GO_input),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  #Save orthogroups IDs to file
  write.table(common_OGs_zscore, paste0(output_dir, "/", my_tissue, "/splsda_zscore_comb/", comb, output_suffix_orthogroups_GO_input),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
