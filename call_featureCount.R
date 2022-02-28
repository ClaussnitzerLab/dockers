library("Rsubread")
args = commandArgs(trailingOnly=TRUE)
input_bam <- args[1]
input_annot_inbuilt <- args[2]
input_gtf <- args[3]
input_featureType <- args[4]
input_attrType <- args[5]
input_strandSpecific <- args[6]
input_isGTFAnnotationFile <- args[7]
input_isPairedEnd <- args[8]

fc <- featureCounts(files=input_bam, annot.inbuilt=input_annot_inbuilt, annot.ext=input_gtf,GTF.featureType=input_featureType, GTF.attrType=input_attrType, strandSpecific=input_strandSpecific, isGTFAnnotationFile=input_isGTFAnnotationFile, isPairedEnd=input_isPairedEnd)
write.csv(fc$counts, "feature_counts.csv")
