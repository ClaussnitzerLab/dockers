library("Rsubread")
args = commandArgs(trailingOnly=TRUE)
input_bam <- args[1]
input_annot_inbuilt <- args[2]
input_gtf <- args[3]
input_featureType <- args[4]
input_attrType <- args[5]
input_strandSpecific <- as.numeric(args[6])
if (args[7] == "TRUE"){
    input_isGTFAnnotationFile <- TRUE
}else{
    input_isGTFAnnotationFile <- FALSE
}

if (args[8] == "TRUE"){
    input_isPairedEnd <- TRUE
}else{
    input_isPairedEnd <- FALSE
}

fc <- featureCounts(files=input_bam, annot.inbuilt=input_annot_inbuilt, annot.ext=input_gtf,GTF.featureType=input_featureType, GTF.attrType=input_attrType, strandSpecific=input_strandSpecific, isGTFAnnotationFile=input_isGTFAnnotationFile, isPairedEnd=input_isPairedEnd)
write.csv(fc$counts, "feature_counts.csv")
