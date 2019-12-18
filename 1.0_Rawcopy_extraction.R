
library(tidyverse)
library(data.table)
library(RODBC)

ch <- odbcConnectAccess2007("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/TP53_DB_v1.accdb") 
# odbcCloseAll()


samples <- sqlFetch(ch, "ELENCO_CAMPIONI_TP53_121119", as.is=TRUE)

SNPS <- samples$CEL_NAME %>% unique %>% na.omit()

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
