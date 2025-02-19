

library(dplyr)
rm(list = ls())

args = commandArgs(trailingOnly=TRUE)
lp_path <- args[1]
variant_path <- args[2]
rsid <- args[3]


root_output <- "LP_output"
dir.create(root_output)

lp <- read.csv(lp_path)
variants <- read.csv(variant_path)
variants <- select(variants, c("ID", rsid))
colnames(variants) <- paste("Metadata_variant", colnames(variants), sep = "_")
df <- merge(lp, variants, by.x = "Metadata_donor", by.y="Metadata_variant_ID", all.x = T)
df_snp <- select(df, c("Metadata_donor", sprintf("Metadata_variant_%s", rsid)))
write.csv(df_snp, sprintf('%s/df_snp.csv', root_output), row.names = F)

compare_alleles <- function(merged, all1, all2){
  # 00 vs 11
  #all1 <- "00"
  #all2 <- "11"
  curr1 <- merged[which(merged$rsID %in% c(all1)),]
  curr2 <- merged[which(merged$rsID %in% c(all2)),]
  
  curr <- merged[which(merged$rsID %in% c(all1, all2)),]
  if (min(curr$LP) < 0){
    curr$LP <- curr$LP-min(curr$LP)
  }
  if (length(unique(curr$rsID)) != 2 | (min(curr$LP) == max(curr$LP)) | nrow(curr1) < 3 | nrow(curr2) < 3){
    tmp <- data.frame(all1, all2, 0, 0, 1, 0,0)
  }else{
    wt <- wilcox.test(curr1$LP, curr2$LP)
    tmp <- data.frame(all1, all2, nrow(curr1), nrow(curr2), wt$p.value,
                      mean(curr[which(curr$rsID == all1), "LP"]), mean(curr[which(curr$rsID == all2), "LP"]))
  }
  colnames(tmp) <- c("Allele1", "Allele2", "N1", "N2", "Wilcox_p", "Mean1", "Mean2")
  return(tmp)
}

# merged <- curr_data
process_gene_var_pair <- function(merged){
  colnames(merged) <- c("IID", "LP", "rsID")
  merged$rsID <- ifelse(merged$rsID == "0/0", "00", 
                        ifelse(merged$rsID == "1/1", "11", "01"))
  # summary(as.factor(merged$rsID))
  
  tmp_00_11 <- compare_alleles(merged, "00", "11")
  tmp_00_01 <- compare_alleles(merged, "00", "01")
  tmp_01_11 <- compare_alleles(merged, "01", "11")
  
  merged_00 <- merged
  merged_00$rsID <- ifelse(merged_00$rsID == "01", "00", merged_00$rsID)
  merged_00$rsID <- ifelse(merged_00$rsID == "00", "0001", merged_00$rsID)
  tmp00_00_11 <- compare_alleles(merged_00, "0001", "11")
  
  merged_11 <- merged
  merged_11$rsID <- ifelse(merged_11$rsID == "01", "11", merged_11$rsID)
  merged_11$rsID <- ifelse(merged_11$rsID == "11", "1101", merged_11$rsID)
  tmp11_00_11 <- compare_alleles(merged_11, "00", "1101")
  
  out <- rbind(tmp_00_11, tmp_00_01, tmp_01_11, tmp00_00_11, tmp11_00_11)
  out <- out[which(out$N1 != 0 & out$N2 != 0),]
  return(out)
}

meta <- select(df, colnames(df)[grepl("Metadata", colnames(df))])
meta_vars <- select(df, colnames(df)[grepl("Metadata_variant", colnames(df))])
lps <- select(df, -colnames(df)[grepl("Metadata", colnames(df))])
a_feat <- colnames(lps)[1]
first <- T
a_var <- colnames(meta_vars)[1]
overall_report <- data.frame(matrix(nrow = 0, ncol = 2))
for (a_var in colnames(meta_vars)){
  for (a_feat in colnames(lps)){
    curr_data <- select(df, c("Metadata_donor", a_feat, a_var))
    
    if (nrow(curr_data) != nrow(curr_data[which(complete.cases(curr_data)),])){
      curr_data <- curr_data[which(complete.cases(curr_data)),]
    }
    
    out <- process_gene_var_pair(curr_data)
    if (nrow(out) > 0){
      out$SNP <- a_var
      out$LP <- a_feat
      out$File <- basename(lp_path)
      if (first){
        first <- F
        overall_report <- out
      }else{
        overall_report <- rbind(out, overall_report)
      }            
    }
  }  
}


write.csv(overall_report, sprintf('%s/overall_report.csv', root_output), row.names = F)




