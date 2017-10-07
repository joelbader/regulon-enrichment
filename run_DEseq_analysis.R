#!/usr/bin/env Rscript
#Set the working directory to the location of this script (necessary for Rstudio):
if(Sys.getenv("RSTUDIO") == "1") setwd(getSrcDirectory(function(x) {x})) 

# Source file containing functions
source("DEseq_analysis_functions.R")

# Get table of Rv to functional categories
Rv_to_category = read.table("all_functional_categories.csv", header=TRUE, sep=",")
#Rv_to_category[,1] = NULL

# Get table of gene name mappings
gene_name_mappings = read.table("V735_to_Rv_MT_annot_mapping.csv", header=FALSE, sep=",", quote="")
gene_name_mappings = merge(gene_name_mappings, Rv_to_category, by.x='V2', by.y='Gene.Name', sort=FALSE, all.x = TRUE)

rownames(gene_name_mappings) = gene_name_mappings[,2]
gene_name_mappings[,2] = NULL
colnames(gene_name_mappings) = c('Rv Homologs (NCBI)','MT Homologs (BLAST)', 'Evalues of MT Homologs (BLAST)','Annotations (NCBI)','Functional.Category')
gene_name_mappings = gene_name_mappings[order(rownames(gene_name_mappings)),]

# Read deseqFile. Expected input is from EDGE-pro
datafile = "deseqFile"
countTable = read.table(datafile, header=TRUE)
countTable = initialize_countTable(countTable) #Use this function to remove all samples that shouldn't factor into this analysis. You could also use it to remove the outliers.
outliers_to_rm = c("DR22","DR25","DR26") # Remove samples with low numbers of mapping reads.
idx_to_rm = c(100000) #Hack to avoid an error due to a NULL when outliers_to_rm is empty (because R sucks)
for(outlier in outliers_to_rm){
  idx_to_rm = c(idx_to_rm, grep(outlier, colnames(countTable), value=FALSE))
}
countTable = countTable[-idx_to_rm]

# Run DESeq2 to define parameters for use in hypothesis testing.
library("DESeq2")
condition = factor( c("7H9", "7H9", "7H9", "PBS", "PBS", "PBS", "lowPhos","lowPhos", "lowPhos", "hypoxia","hypoxia","rpoB-7H9", "rpoB-7H9", "rpoB-7H9", "rpoB-7H9", "rpoB-PBS", "rpoB-PBS", "rpoB-lowPhos", "rpoB-lowPhos", "rpoB-lowPhos", "rpoB-hypoxia","rpoB-hypoxia"))
condition = condition[-idx_to_rm]
colData = data.frame(condition)
row.names(colData) = colnames(countTable)
dds = DESeqDataSetFromMatrix(countData = countTable, colData = colData, design = ~ condition)
dds = DESeq(dds)

#Run analysis pipeline
list[res7H9vshyp, res7H9vsPBS, res7H9vslowPhos] = wj_genotype_fixed("")

## Output analysis to file
#_low-read-rm
write.csv(res7H9vshyp, file = "Analysis_Output/7H9vshyp_low-read-rm.csv")
write.csv(res7H9vsPBS, file = "Analysis_Output/7H9vsPBS_low-read-rm.csv")
