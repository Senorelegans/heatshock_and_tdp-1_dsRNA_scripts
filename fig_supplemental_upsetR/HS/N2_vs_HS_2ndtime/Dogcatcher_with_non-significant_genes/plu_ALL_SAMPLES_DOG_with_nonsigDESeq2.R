

##########################################################################################
##########################################################################################
########DESeq2

library(DESeq2)
library("BiocParallel")
register(SerialParam())
register(MulticoreParam(4))

##Read in the count matrix
cts <- as.matrix(read.csv("HEATSHOCK_LFC_FLIPPED_DOG_2ndtime/N2_vs_HS_2ndtime/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_Rsubread_sense_normalized_rRNA_star.txt",sep="\t",row.names=1))
coldata <- read.table("HEATSHOCK_LFC_FLIPPED_DOG_2ndtime/initial_Rsubread/col_data.txt")

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
write.csv(res, file="HEATSHOCK_LFC_FLIPPED_DOG_2ndtime/N2_vs_HS_2ndtime/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_DESeq2_sense.csv")
#############################################





cts <- as.matrix(read.csv("HEATSHOCK_LFC_FLIPPED_DOG_2ndtime/N2_vs_HS_2ndtime/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_Rsubread_antisense_normalized_rRNA_star.txt",sep="\t",row.names=1))



#LRT test
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ assay + condition + assay:condition)


sizeFactors(dds)=mylist$SF    #Makes the size factors for dds from the col_data length


dds <- dds[ rowSums(counts(dds)) > 2,]   #Take everything over 2 because we added one in the normalization python steps

dds$condition <- factor(dds$condition, levels=c("T","C"))
dds <- DESeq(dds, test="LRT", reduced= ~ assay + condition, parallel = TRUE)
res <- results(dds)
write.csv(res, file="HEATSHOCK_LFC_FLIPPED_DOG_2ndtime/N2_vs_HS_2ndtime/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_DESeq2_antisense.csv")
#############################################
