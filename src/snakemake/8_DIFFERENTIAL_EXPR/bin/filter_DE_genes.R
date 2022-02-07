#### This is a script to compute GO enrichments with a custom annotation

#upload libraries
library("tidyverse")
library("reshape2")

###### Read arguments
Args = commandArgs(trailingOnly=TRUE);
expr_df_file = (Args[1])
DE_genes_file = (Args[2])
my_tissue = (Args[3])
vertebrata_species = (Args[4])
insecta_species = (Args[5])
outgroup_species = (Args[6])
log2fc_cutoff = (Args[7])
pvalue_cutoff = (Args[8])
output_file = (Args[9])

###### Define variables
Vertebrata = strsplit(vertebrata_species, ",")[[1]] 
Insecta = strsplit(insecta_species, ",")[[1]]
Outgroups = strsplit(outgroup_species, ",")[[1]]
log2fc_cutoff = as.numeric(log2fc_cutoff)
pvalue_cutoff = as.numeric(pvalue_cutoff)


###### Load the expr df at the tissue level.
expr_df = read.table(expr_df_file, header=TRUE, row=1)

###### Load dataframe with differentially expressed genes
diff_expr_OGs_df = read.delim(DE_genes_file, header=TRUE, row=1)
diff_expr_OGs_df = subset(diff_expr_OGs_df, log2FoldChange >= log2fc_cutoff & padj <= pvalue_cutoff)
diff_expr_OGs = rownames(diff_expr_OGs_df)

###### Compute zscores among the median expression.
diff_expr_df = expr_df[rownames(expr_df) %in% diff_expr_OGs,]
diff_expr_df_long = melt(as.matrix(diff_expr_df))
colnames(diff_expr_df_long) = c("OG_ID", "Metasample", "Expr")
diff_expr_df_long$Species = sub("_.*", "", diff_expr_df_long$Metasample)
diff_expr_df_long$Tissue = sub(".*_", "", diff_expr_df_long$Metasample)
#Add clade so that the dataframe could be grouped by that.
clade_expr_df = diff_expr_df_long %>% 
  mutate(Clade = case_when(Species %in% Vertebrata ~ "Vertebrata",
                           Species %in% Insecta ~ "Insecta",
                           Species %in% Outgroups ~ "Outgroups")
               )
      
OG_zscore_df = clade_expr_df %>% group_by(OG_ID, Tissue, Clade) %>% 
  summarize(median_expr=median(Expr, na.rm=TRUE)) %>%
  ungroup %>%
  group_by(OG_ID) %>% 
  mutate(zscore = (median_expr - mean(median_expr, na.rm=TRUE))/sd(median_expr, na.rm=TRUE))

####### Require the zscore in the analyzed tissue to be higher than 1 in at least Vertebrate and Insects
#OG_to_select_counts = table(subset(OG_zscore_df, Clade %in% c("Vertebrata", "Insecta") & Tissue==my_tissue & zscore >= 1)$OG_ID)
#OG_to_select = names(OG_to_select_counts[OG_to_select_counts == 2])
OG_to_select_counts = table(subset(OG_zscore_df, Clade %in% c("Vertebrata", "Insecta", "Outgroups") & Tissue==my_tissue & zscore >= 1)$OG_ID)
OG_to_select = names(OG_to_select_counts[OG_to_select_counts == 3])


####### Filter file with the differentially expressed genes
filtered_diff_expr_OGs = diff_expr_OGs_df[OG_to_select,]
####### Save to file
write.table(filtered_diff_expr_OGs, output_file, sep="\t", quote=FALSE, col.names=NA)
