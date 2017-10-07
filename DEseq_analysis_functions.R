# List trick to allow for multiple returned values
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
   args <- as.list(match.call())
   args <- args[-c(1:2,length(args))]
   length(value) <- length(args)
   for(i in seq(along=args)) {
      a <- args[[i]]
      if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
   }
   x
}

#Identify genes found to be significantly up and down regulated across all stress conditions relative to 7H9 (rich media)
wj_genotype_fixed<-function(genotype)
{
   test_thresh = 0

   seven_H9 = paste(genotype, "7H9", sep="")
   hypoxia = paste(genotype, "hypoxia", sep="")
   PBS = paste(genotype, "PBS", sep="")
   lowPhos = paste(genotype, "lowPhos", sep="")

   r7H9vshyp = results(dds, contrast=c("condition", hypoxia, seven_H9), lfcThreshold = test_thresh, altHypothesis = "greaterAbs")
   r7H9vsPBS = results(dds, contrast=c("condition", PBS, seven_H9), lfcThreshold = test_thresh, altHypothesis = "greaterAbs")
   r7H9vslowPhos = results(dds, contrast=c("condition", lowPhos, seven_H9), lfcThreshold = test_thresh, altHypothesis = "greaterAbs")

   direction_hyp = r7H9vshyp[,2] # Initialize to the proper length
   direction_hyp[r7H9vshyp[,2] > 0] = "UP"
   direction_hyp[r7H9vshyp[,2] < 0] = "DOWN"
   r7H9vshyp = transform(merge(gene_name_mappings, data.frame(r7H9vshyp[,c(1,2,3,5,6)],direction_hyp), by='row.names', sort=FALSE), row.names=Row.names, Row.names = NULL)
   colnames(r7H9vshyp)[colnames(r7H9vshyp) == 'direction_hyp'] = 'Expression Direction (Hypoxia/7H9)'

   direction_pbs = r7H9vsPBS[,2] # Initialize to the proper length
   direction_pbs[r7H9vsPBS[,2] > 0] = "UP"
   direction_pbs[r7H9vsPBS[,2] < 0] = "DOWN"
   r7H9vsPBS = transform(merge(gene_name_mappings, data.frame(r7H9vsPBS[,c(1,2,3,5,6)],direction_pbs), by='row.names', sort=FALSE), row.names=Row.names, Row.names = NULL)
   colnames(r7H9vsPBS)[colnames(r7H9vsPBS) == 'direction_pbs'] = 'Expression Direction (PBS/7H9)'

   direction_lowPhos = r7H9vslowPhos[,2] # Initialize to the proper length
   direction_lowPhos[r7H9vslowPhos[,2] > 0] = "UP"
   direction_lowPhos[r7H9vslowPhos[,2] < 0] = "DOWN"
   r7H9vslowPhos = transform(merge(gene_name_mappings, data.frame(r7H9vslowPhos[,c(1,2,3,5,6)],direction_lowPhos), by='row.names', sort=FALSE), row.names=Row.names, Row.names = NULL)
   colnames(r7H9vslowPhos)[colnames(r7H9vslowPhos) == 'direction_lowPhos'] = 'Expression Direction (Low Phosphate/7H9)'

   return(list(r7H9vshyp, r7H9vsPBS, r7H9vslowPhos))
}

#Initializes countTable so that format and names are correct
initialize_countTable<-function(countTable)
{
   # Rename any duplicate gene names. Add _1, _2, ...
   gene_list = as.character(countTable[,1])
   dup_index = which(duplicated(gene_list))
   for (i in 1:length(dup_index)){
      dup_genename = gene_list[dup_index[i]]
      n = length(which(gene_list == dup_genename))
      gene_list[gene_list == dup_genename] = paste(gene_list[gene_list == dup_genename], '_', 1:n, sep='')
   }
   # Remove first row of dataframe and use row.names instead
   countTable[,1] = as.factor(gene_list)
   rownames(countTable) = countTable[,1]
   countTable[,1] = NULL

   # Remove the EmptyGene and EmptyRNA included on each supercontig (necessary to avoid bug in EDGE-pro)
   countTable = countTable[!rownames(countTable) %in% c("EmptyRNA_1","EmptyGene_1","EmptyRNA_2","EmptyGene_2"), ]

   return(countTable)
}
