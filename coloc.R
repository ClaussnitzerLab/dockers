library(arrow)
library(coloc)
library(stringr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
cell_state <- args[1]
lp_root <- args[2]
lp_feat_path <- args[3]
trait_path <- args[4]
chr_index_input <- args[5]
gwas_n <- args[6]


lp_features <- read_parquet(lp_feat_path)
trait <- read.csv(trait_path);

tmp <- as.data.frame(t(as.data.frame(str_split(trait$chrbpa1a2_hg38, ":"))))
trait$coord <- paste(tmp$V1, tmp$V2, sep = ":")
trait <- select(trait, c("coord", "P"))
trait <- trait %>% group_by(coord) %>% summarise(P = min(P))
colnames(trait) <- c("coord", "trait_p")

chr_index <- 22
first <- T
for (chr_index in c(chr_index_input)){
  fpath <- sprintf("%s/%s_chr%d.parquet", lp_root, cell_state, chr_index)
  lp <- read_parquet(fpath); head(lp) 
  curr_trait <- trait
  for (a_phen in unique(lp_features$lps)){
    
    curr_gwas <- lp[which(lp$phenotype_id == a_phen),]
    tmp <- as.data.frame(t(as.data.frame(str_split(curr_gwas$variant_id, ":"))))
    curr_gwas$coord <- paste(tmp$V1, tmp$V2, sep = ":")
    curr_gwas$af <- ifelse(curr_gwas$af < .5, curr_gwas$af, 1-curr_gwas$af)
    curr_gwas <- select(curr_gwas, c("coord", "pval", "af"))
    colnames(curr_gwas) <- c("coord", "lp_p", "af")
    
    input_var_id <- merge(curr_trait, curr_gwas, by="coord")
    if (nrow(input_var_id) > 1){
      result <- coloc.abf(dataset1=list(pvalues=input_var_id$trait_p, snp=input_var_id$coord, type="quant", N=gwas_n),
                  dataset2=list(pvalues=input_var_id$lp_p, snp=input_var_id$coord, type="quant", N=118), 
                  MAF=input_var_id$af)
      tmp <- data.frame(a_phen, chr_index, result$summary[1], result$summary[2], result$summary[3],
                        result$summary[4], result$summary[5], result$summary[6])
      colnames(tmp) <- c("phenotype", "chr_index", "nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
                         
    tmp_snp <- as.data.frame(result$results)
    tmp_snp$phenotype <- a_phen
    tmp_snp$chr_index <- chr_index
    tmp_snp <- tmp_snp[which(tmp_snp$SNP.PP.H4 > .01),]


      if(first){
        first <- F
        out <- tmp
        out_snp <- tmp_snp
      }else{
        out <- rbind(out, tmp)
        out_snp <- rbind(out_snp, tmp_snp) 
      }
    }
  }
}    

writeLines(sprintf("Num row: %d", nrow(out)))
write_parquet(out, "out.parquet", row.names=F)
write_parquet(out_snp, "out_snp.parquet", row.names=F)

