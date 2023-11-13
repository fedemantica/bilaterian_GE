#This is to compute the zscore by species
######## Upload libraries
library("reshape2")
library("tidyverse")

######## Set Arguments
Args <- commandArgs(trailingOnly=TRUE);
expr_table = (Args[1]);
output_file = (Args[2]);

######## Main
#Read inputs
all_OG_df = read.delim(expr_table, header=TRUE, row=1)
#Transform to long format
all_OG_long_df = melt(as.matrix(all_OG_df))
colnames(all_OG_long_df) = c("OG_ID", "Metasample", "Expr")
all_OG_long_df$Species = sub("_.*", "", all_OG_long_df$Metasample)

#Compute Zscores by species for each conserved genes
zscore_long_df = all_OG_long_df %>% 
  group_by(Species, OG_ID) %>%
  mutate(Zscore = (Expr - mean(Expr, na.rm=TRUE))/sd(Expr, na.rm=TRUE)) %>%
  ungroup() %>%
  select(-c(Expr, Species))

#Transform to wide format
zscore_df = zscore_long_df %>%
  pivot_wider(names_from=Metasample, values_from=Zscore)

zscore_df = as.data.frame(zscore_df)
rownames(zscore_df) = zscore_df$OG_ID
zscore_df$OG_ID = NULL

#Save to output
write.table(zscore_df, output_file, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)