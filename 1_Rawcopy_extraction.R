
library(tidyverse)
library(data.table)

samples <- fread("D:/analisi_in_corso/TP53/ELENCO_CAMPIONI_TP53_definitivo_220719.txt") 

SNPS <- samples$CEL_NAME %>% unique

rawcopydir <- "D:/RAW_DATA/RAWCOPY/alpha=10^-7/"


results <- data.frame()

for(i in SNPS){
  
  try({
  print(i)
  
  setwd(paste0(rawcopydir,"/",i))
  
  segms <- data.table::fread("segments2.txt") 
  names(segms)[5] <- "logR"
  names(segms)[7] <- "Allelic_imbalance"
  
  namedsegms <- cbind(sample=i, segms)
  
  results <- rbind(results,namedsegms)
  })
}

readr::write_tsv(results, "D:/analisi_in_corso/TP53/CN_calls/rawcopy_segm_TP53_cohort.txt")
