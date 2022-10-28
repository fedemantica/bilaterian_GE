#### This is a script to compute GO enrichments with a custom annotation

#upload libraries
library("gprofiler2")

#read arguments
Args = commandArgs(trailingOnly=TRUE);
gene_set_file = (Args[1]) 
background_file = (Args[2])
annotation_file = (Args[3])
output_file = (Args[4])

#gene_set_file = "/Users/federica/mnt/projects/ludo_analysis/data/GO_enrichments/gene_sets/microexons/hg38_gene_set.txt"
#background_file = "/Users/federica/mnt/projects/ludo_analysis/data/GO_enrichments/GO_background/hg38_background.txt"
#background_file = "/Users/federica/mnt/projects/ludo_analysis/data/GO_backgrounds/Hs2_RetinaExpr_average_5cRPKM_NORM.txt"
#backgroundfile = "/Users/federica/mnt/projects/ludo_analysis/data/GO_enrichments/GO_background/hg38_background.txt"
#annotation_file = "/Users/federica/mnt/projects/ludo_analysis/data/GO_enrichments/GO_annot/test.gmt"
#output_file = "/Users/federica/Documents/CRG_PhD/weekly_meetings/21_07_22/test_GO.tab"

my_genes = as.vector(read.delim(gene_set_file, header=FALSE, col.names=c("genes"))$genes)
my_background_genes = as.vector(read.delim(background_file, header=FALSE, col.names=c("genes"))$genes)
token = upload_GMT_file(gmtfile = annotation_file)
#custom_gp = gost(my_genes, organism = token, correction_method="fdr", custom_bg = my_background_genes, significant=TRUE)
custom_gp = gost(my_genes, organism = token, correction_method="fdr", custom_bg = my_background_genes, significant=FALSE)
my_res = custom_gp$result
#save to output file
write.table(my_res[,1:12], output_file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
