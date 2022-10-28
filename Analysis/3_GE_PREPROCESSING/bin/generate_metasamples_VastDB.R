#This is a script to derive 3 groups of samples within each tissue
#In principle it can be applied only when more than 3 samples are actually present

#It requires as input the sva-log2-TPM corrected tables and the relative metatable
#NB: I am using the MEDIAN right now

#Define arguments
Args = commandArgs(trailingOnly=TRUE);
expr_file = (Args[1]) #sva corrected tables
metadata_file = (Args[2])
expr_output = (Args[3]) #file with average expressions across groups of samples
metadata_output = (Args[4]) #file with the correspondence samples - group. The groups will be named Tissue_num (ex: Adipose_1, Adipose_2, Adipose_3)

#For debugging
# species = "Spu"
# expr_file = paste0("/Users/federica/mnt/projects/bilaterian_GE/data/preprocessing/sva_correction/", species, "-selected_samples_gene_SVA-log2-TPMs.tab")
# metadata_file = paste0("/Users/federica/mnt/projects/bilaterian_GE/data/samples_metadata/", species, "_samples_info.tab")
# expr_output = paste0("/Users/federica/mnt/projects/bilaterian_GE/data/preprocessing/metasamples/", species, "-metasamples_median_expr.tab")
# metadata_output = paste0("/Users/federica/mnt/projects/bilaterian_GE/data/preprocessing/metasamples/", species, "-metasamples_metadata.txt")

#Read inputs
expr_table = read.delim(expr_file, header=TRUE, row=1)
metadata_table = read.delim(metadata_file, header=TRUE, row=4)

tissues = unique(as.vector(metadata_table$Tissue)) #Define tissues on which to cycle based on metatable
#Cycle on the tissues
all_metasample_expr_table = as.data.frame(rownames(expr_table)); colnames(all_metasample_expr_table) = "GeneID" #Ugly trick to make the cbind work
######
#Test for Cdi and Bge
rownames(all_metasample_expr_table) = rownames(expr_table)
######
all_metasample_metadata_table = vector()

for (tissue in tissues) {
  tissue_samples = rownames(subset(metadata_table, Tissue==tissue))
  tissue_samples_selected = tissue_samples[tissue_samples %in% colnames(expr_table)] #Get only the samples that were not removed by the sva iteration correction
  tissue_expr_table = expr_table[,tissue_samples_selected] #Subset expression table to only tissue samples
  
  #cycle on the VastDB tissue groups
  tissue_groups = unique(as.vector(subset(metadata_table, Tissue==tissue)$Group)) #Isolate VastDB metasamples for each tissue
  for (metasample in tissue_groups) {
    metasample_samples = rownames(subset(metadata_table, Group==metasample))
    metasample_samples = metasample_samples[metasample_samples %in% tissue_samples_selected]
    if (length(metasample_samples) == 0) {
      next
      } else if (length(metasample_samples) >=2) {
        metasample_expr_table = tissue_expr_table[,metasample_samples]
        metasample_avg_expr_table = as.data.frame(apply(metasample_expr_table, 1, median)) #Compute gene median between all samples in the metasample
        colnames(metasample_avg_expr_table) = metasample #Add tissue_metasample as colname
        all_metasample_expr_table = cbind(all_metasample_expr_table, 
                                          metasample_avg_expr_table) #join to general median expr table
        metasample_metadata_table = data.frame(Sample = metasample_samples, 
                                               Metasample = rep(metasample, length(metasample_samples)),
                                               Tissue = rep(tissue, length(metasample_samples))) #generate metadata table
        all_metasample_metadata_table = rbind(all_metasample_metadata_table, metasample_metadata_table)
      } else { #If there is only one sample per tissue group
        if (!is.data.frame(tissue_expr_table)) { #if there is only one sample PER tissue, I need to retransform into a dataframe
          metasample_expr_table = as.data.frame(tissue_expr_table)
        } else {
          metasample_expr_table = as.data.frame(tissue_expr_table[,metasample_samples])
        }
        colnames(metasample_expr_table) = metasample
        rownames(metasample_expr_table) =  all_metasample_expr_table$GeneID #Fede added 27/08/21
        all_metasample_expr_table = cbind(all_metasample_expr_table, 
                                          metasample_expr_table) #Join to general median expr table
        metasample_metadata_table = data.frame(Sample = metasample_samples, 
                                               Metasample = metasample,
                                               Tissue = tissue)
        all_metasample_metadata_table = rbind(all_metasample_metadata_table, metasample_metadata_table)
    } 
  }
}
  
all_metasample_expr_table$GeneID = NULL #Ugly trick to make the cbind work

#Save to output files
write.table(all_metasample_expr_table, expr_output, sep="\t", quote=FALSE, row.names = TRUE, col.names=NA)
write.table(all_metasample_metadata_table, metadata_output, sep="\t", quote=FALSE, row.names = FALSE, col.names=TRUE)