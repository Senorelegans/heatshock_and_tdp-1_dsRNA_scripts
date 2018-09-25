#########################################################################################
#####Run Rsubread for counts
library(Rsubread);

files=c(bamlist)

gtf=("initial_annotation_file")

###Sense Unique
sense=featureCounts(files,
isGTFAnnotationFile = TRUE,
annot.ext = gtf,
GTF.attrType = "gene_id",
allowMultiOverlap = TRUE,
nthreads = cpus,
countMultiMappingReads=TRUE,
strandSpecific = 1)

write.table(x=data.frame(
sense$annotation[,c("GeneID","Length")],
sense$counts,stringsAsFactors=FALSE),
file="outputprefix/Rsubread_sense.txt",
quote=FALSE,sep="\t",row.names=FALSE)

###Anti-sense 
antisense=featureCounts(files,
isGTFAnnotationFile = TRUE,
annot.ext = gtf,
GTF.attrType = "gene_id",
allowMultiOverlap = TRUE,
countMultiMappingReads=TRUE,
nthreads = cpus,
strandSpecific = 2)

write.table(x=data.frame(
antisense$annotation[,c("GeneID","Length")],
antisense$counts,stringsAsFactors=FALSE),
file="outputprefix/Rsubread_antisense.txt",
quote=FALSE,sep="\t",row.names=FALSE)
