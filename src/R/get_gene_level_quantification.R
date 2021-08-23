#Rscript to convert transcript-wise to gene-wise gene expression quantification
#saving both TPMs and counts

#Upload libraries
library("tximportData")
library("tximport")

#set arguments
Args = commandArgs(trailingOnly=TRUE);
tx2gene_file = (Args[1]) 
output_TPMs = (Args[2])
output_counts = (Args[3])

files = vector()
for (index in seq(4, length(Args))) {
  files = c(files, (Args[index]))
}
#SE_downloaded_files = (Args[4])
#SE_in_house_files = (Args[5]) 
#PE_downloaded_files = (Args[6]) 
#PE_in_house_files = (Args[7])

#building file list
#files = c(SE_downloaded_files, SE_in_house_files, PE_downloaded_files, PE_in_house_files)
names(files) = basename(sub("abundance.tsv", "", files))

#debugging
#print(class(SE_downloaded_files)); print(length(SE_downloaded_files))
print(files)
print(length(files))

#uploading tx2gene file
tx2gene = read.delim(tx2gene_file, header=TRUE)
#operating the conversion
txi.kallisto.tsv = tximport(files, type="kallisto", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
#isolating tables with gene-level quantification
gene_counts_table = txi.kallisto.tsv$counts
gene_TPMs_table = txi.kallisto.tsv$abundance
print(txi.kallisto.tsv$countsFromAbundance)
print(attributes(txi.kallisto.tsv))
# print(gene_counts_table)
# print(gene_TPMs_table)
#print(txi.kallisto.tsv$abundance)
#print(sum(txi.kallisto.tsv$countsFromAbundance[,2]))

#saving to output file
write.table(gene_counts_table, output_counts, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
write.table(gene_TPMs_table, output_TPMs, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
