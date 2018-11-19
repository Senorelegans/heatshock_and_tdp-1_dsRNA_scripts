#########################################################################################
#####Run Rsubread for counts
library(Rsubread);

files=c('BAMS/N2WT_INP_.bam', 'BAMS/OK803-INP__.bam', 'BAMS/N2-J2-1_.bam', 'BAMS/N2-J2-2_.bam', 'BAMS/N2-J2-3_.bam', 'BAMS/OK803-J2-1_.bam', 'BAMS/OK803-J2-2_.bam', 'BAMS/OK803-J2-3_.bam')

gtf=("OK_LFC_FLIPPED_DOG_2ndtime/N2_vs_OK_2ndtime/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_DOG_with_nonsig.gtf")

###Sense Unique
sense=featureCounts(files,
isGTFAnnotationFile = TRUE,
annot.ext = gtf,
GTF.attrType = "gene_id",
allowMultiOverlap = TRUE,
nthreads = 4,
countMultiMappingReads=TRUE,
strandSpecific = 1)

write.table(x=data.frame(
sense$annotation[,c("GeneID","Length")],
sense$counts,stringsAsFactors=FALSE),
file="OK_LFC_FLIPPED_DOG_2ndtime/N2_vs_OK_2ndtime/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_Rsubread_sense.txt",
quote=FALSE,sep="\t",row.names=FALSE)

###Anti-sense 
antisense=featureCounts(files,
isGTFAnnotationFile = TRUE,
annot.ext = gtf,
GTF.attrType = "gene_id",
allowMultiOverlap = TRUE,
countMultiMappingReads=TRUE,
nthreads = 4,
strandSpecific = 2)

write.table(x=data.frame(
antisense$annotation[,c("GeneID","Length")],
antisense$counts,stringsAsFactors=FALSE),
file="OK_LFC_FLIPPED_DOG_2ndtime/N2_vs_OK_2ndtime/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_Rsubread_antisense.txt",
quote=FALSE,sep="\t",row.names=FALSE)
