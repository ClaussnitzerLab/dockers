#!/usr/bin/env Rscript
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

write.csv(as.data.frame(args), "feature_counts.csv")
