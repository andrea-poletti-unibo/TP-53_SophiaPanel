library(data.table)
library(tidyverse)

#=========== merge FAILED arrays results ============
failed1 <- fread("../ASCAT/BATCH1_6.0/BATCH1_6.0_failed_arrays_list.csv")
failed2 <- fread("../ASCAT/BATCH2_HD/BATCH2_HD_failed_arrays_list.csv")
failed3 <- fread("../ASCAT/BATCH3_HD/BATCH3_HD_failed_arrays_list.csv")

failed_all <- rbind(failed1[-1,-1], failed2[-1,-1], failed3[-1,-1])
failed <- failed_all$V2 %>% str_replace("\\.GenomeWideSNP_6\\.","(GenomeWideSNP_6)") %>% str_replace("\\.CytoScanHD_Array\\.","(CytoScanHD_Array)") %>% str_replace("^X","")
failed <- data.frame(sample = failed)

write_tsv(failed, "../ASCAT/All_failed_samples.txt")

#========== merge batch SEGMENTS results =============

segm1 <- fread("../ASCAT/BATCH1_6.0/BATCH1_6.0adj_fitted_and_raw_segments.tsv")
segm2 <- fread("../ASCAT/BATCH2_HD/BATCH2_HDadj_fitted_and_raw_segments.tsv")
segm3 <- fread("../ASCAT/BATCH3_HD/BATCH3_HDadj_fitted_and_raw_segments.tsv")

segm_all <- rbind(segm1,segm2,segm3)

segm_all$sample <- segm_all$sample %>% str_replace("\\.GenomeWideSNP_6\\.","(GenomeWideSNP_6)") %>% str_replace("\\.CytoScanHD_Array\\.","(CytoScanHD_Array)") %>% str_replace("^X","")

write_tsv(segm_all, "../ASCAT/All_adj_fitted_and_raw_segments.tsv")

