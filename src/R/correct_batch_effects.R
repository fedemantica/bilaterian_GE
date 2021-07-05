#Rscript to be run from a snakemake
#Here I aim at removing samples outlier and correct hidden batch effects
#I am referring to the procedure developed by Fukushima et al, NatCom, 2020.
#I slightly modify it so that the samples from the same Bioproject can be independently removed one from the other.

#upload libraries
library("sva")
library("limma")
library("tidyverse")
library("hashmap")

#define variables from each argument
Args = commandArgs(trailingOnly=TRUE);
log2_TPM_file = (Args[1]) 
metadata_file = (Args[2]) 
output_file = (Args[3]) #path to output file

#define function to correct the expression values based on the identified surrogate variables.
#function from: https://support.bioconductor.org/p/47350/
cleaningY = function(y, mod, svaobj) {
  X=cbind(mod,svaobj$sv)
  Hat=solve(t(X)%*%X)%*%t(X)
  beta=(Hat%*%t(y))
  P=ncol(mod)
  cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}

#read input dataframe
#input_df = read.delim(log2_TPM_file, header=TRUE, sep="\t")
input_df_raw = read.delim("/Users/federica/mnt/projects/bilaterian_GE/data/preprocessing/kallisto_out/Xtr/all_samples_TPMs.tab", 
                      header=TRUE, sep="\t", row=1)
input_df_raw = select(input_df_raw, -c("Adipose_a", "Muscle_a"))
#I NEED TO UPLOAD LOG2 DATA.
#I AM TEMPORARILY COMPUTING THE LOG2 HERE
input_df = apply(input_df_raw, c(1,2), function(x) log2(x+1))

#transform dataframe into a matrix
input_matrix = as.matrix(input_df)
#removing the genes where expression is 0 in all the samples
row_sub = apply(input_matrix, 1,function(x) !all(x==0))
input_matrix_noZeros = input_matrix[row_sub,]

#read input metadata
#metadata_df = read.delim(metadata_file, header=TRUE, sep="\t", row=4)
metadata_df = read.delim("/Users/federica/mnt/projects/bilaterian_GE/data/samples_metadata/Xtr_samples_info.tab", 
                         header=TRUE, sep="\t", row=4)


metadata_df = subset(metadata_df, !(Tissue %in% c("Adipose", "Muscle")))
metadata_df$Tissue = droplevels(metadata_df$Tissue)
#build full model necessary for the sva function
full_model = model.matrix(~as.factor(Tissue), metadata_df)
#build null model necessary for the sva function (in this case only intercept: we are using surrogate variables identified by the model)
null_model = model.matrix(~1, metadata_df)
                    
#apply sva function. Let the model estimate the number of necessary surrogate variables.
svobj = sva(input_matrix_noZeros, full_model, null_model)
            
#correct the batch effects with the removeBatchEffect function from Limma
corrected_expr_matrix = cleaningY(input_matrix_noZeros, full_model, svobj)

#calculate the average expression among all samples of the same tissue group (one value per gene)
corrected_expr_tb = as_tibble(corrected_expr_matrix) %>% mutate(genes=rownames(corrected_expr_matrix))
corrected_expr_tb_long = pivot_longer(corrected_expr_tb, cols=colnames(corrected_expr_matrix))
#add the relative tissue
sample_tissue_dict = hashmap(rownames(metadata_df), as.character(metadata_df$Tissue))
corrected_expr_tb_long$Tissue = sample_tissue_dict$find(corrected_expr_tb_long$name)
#compute the tissue average expression
avg_expr_tb_long = group_by(corrected_expr_tb_long, genes, Tissue) %>% summarise(avg_expr=mean(value))
#Transform to wide format
avg_expr_tb = pivot_wider(avg_expr_tb_long, names_from = Tissue, values_from = avg_expr)
#Add the averages across tissues to the samples table to compute the Pearson correlation
joint_expr_tb = inner_join(x=corrected_expr_tb,
           y=avg_expr_tb,
           by="genes")


#calculate the Pearson correlation between the expression in each samples and in each of the averaged tissue groups
samples_to_remove = vector()
for (sample in rownames(metadata_df)) {
  sample_expr = as.vector(unlist(joint_expr_tb[,sample]))
  sample_tissue = sample_tissue_dict$find(sample) #tissue from which the sample was taken
  tissue_pearson_cor_vector = vector()
  for (tissue in unique(as.character(metadata_df$Tissue))) {
    tissue_expr = as.vector(unlist(joint_expr_tb[,tissue]))
    pearson_cor = cor.test(sample_expr, tissue_expr, method="pearson")
    tissue_pearson_cor_vector = c(tissue_pearson_cor_vector, pearson_cor$estimate)
  }
  names(tissue_pearson_cor_vector) = unique(as.character(metadata_df$Tissue))
  if (max(tissue_pearson_cor_vector) != tissue_pearson_cor_vector[sample_tissue]) { #if the max correlation is not in the tissue from which the samples was originated
    samples_to_remove = c(samples_to_remove, sample)
  }
}

#remove the samples where the maximum Pearson correlation is not detected in the tissue of origin
if (lenght(samples_to_remove)) == 0 {break}
else {}

#repeat the operation until no samples outliers are detected.
#I think this implies starting again from the uncorrected values.





