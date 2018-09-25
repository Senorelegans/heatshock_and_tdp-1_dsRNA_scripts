

##########################################################################################
##########################################################################################
########DESeq2

library(DESeq2)
library("BiocParallel")
register(SerialParam())
register(MulticoreParam(cpus))

##Read in the count matrix
cts <- as.matrix(read.csv("outputprefix/Rsubread_sense_normalized_rRNA_star.txt",sep="\t",row.names=1))
coldata <- read.table("outputprefix/COL_DATA")

#########################################
#      Make your Count data match col Input
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))
#########################################

##########################################
#Initialize list of 1's for size factors
coldata_rows = nrow(coldata)                      #Get number of rows
mylist = list (SF = integer(coldata_rows) )       #Make list of 0's
mylist$SF = mylist$SF + 1                         #Add one to list
##########################################




#LRT test
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ assay + condition + assay:condition)


sizeFactors(dds)=mylist$SF    #Makes the size factors for dds from the col_data length


dds <- dds[ rowSums(counts(dds)) > 2,]   #Take everything over 2 because we added one in the normalization python steps

dds$condition <- factor(dds$condition, levels=c("T","C"))
dds <- DESeq(dds, test="LRT", reduced= ~ assay + condition, parallel = TRUE)
res <- results(dds)
write.csv(res, file="outputprefix/DESeq2_sense.csv")
#############################################





cts <- as.matrix(read.csv("outputprefix/Rsubread_antisense_normalized_rRNA_star.txt",sep="\t",row.names=1))



#LRT test
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ assay + condition + assay:condition)


sizeFactors(dds)=mylist$SF    #Makes the size factors for dds from the col_data length


dds <- dds[ rowSums(counts(dds)) > 2,]   #Take everything over 2 because we added one in the normalization python steps

dds$condition <- factor(dds$condition, levels=c("T","C"))
dds <- DESeq(dds, test="LRT", reduced= ~ assay + condition, parallel = TRUE)
res <- results(dds)
write.csv(res, file="outputprefix/DESeq2_antisense.csv")
#############################################
