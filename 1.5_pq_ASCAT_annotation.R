library(data.table)
library(tidyverse)
library(GenomicRanges)

segms <- fread("../ASCAT/All_segms_refit_center.txt") # import the ASCAT segments

# reformat
segms_prep <- segms %>% select(1,2, start=startpos, end=endpos, CN=CNvalueRaw_WGDref_cen) 

chrdata <- fread("D:/Cartelle_personali/Andrea/GIT_Shirke019/Utility/chrom_arm_data_v2.txt") # import the chromosome arms data (previously generated file)

# transform in GRanges
chrdataGR <- makeGRangesFromDataFrame(chrdata, keep.extra.columns = T)

# create a list of samples names
SAMPLES <- unique(segms$sample)

#____________ Loop for each SAMPLE _____________
cohort_res <- data.frame()
i=1
for (i in 1:length(SAMPLES)){
  print(i)
  print(SAMPLES[i])
  
  segms1 <- filter(segms_prep, sample== SAMPLES[i])
  segms1GR <- makeGRangesFromDataFrame(segms1, keep.extra.columns = T)
  
  #_______ Loop for each CHR ARM ___________
  samp_res <- data.frame()
  j=1
  for (j in 1:nrow(chrdata)) {
    arm <- chrdataGR[j]$chrarm
    
    armInt <- pintersect(segms1GR, chrdataGR[j]) %>% as.data.frame() %>% filter(hit==TRUE)
    
    armRes <- cbind(arm, armInt[,-c(10,12)]) %>% dplyr::rename(chr = seqnames) %>% select(7,2:5,1,8)
    
    samp_res <- rbind(samp_res, armRes)
  }
  
  cohort_res <- rbind(cohort_res, samp_res)
}

write_tsv(cohort_res, "../ASCAT/All_segmsRC_pq.txt")



#================= IGV files ================

IGV_segms_prep <- segms_prep %>% mutate(num.mark=NA) %>% select(ID=sample, chrom=chr, loc.start=start, loc.end= end, num.mark, seg.mean=CN)
write_tsv(IGV_segms_prep, "../ASCAT/IGV_files/IGV_segmsRC_TP53.seg")

IGV_cohort_res <- cohort_res %>% mutate(num.mark=NA) %>% select(ID=sample, chrom=chr, loc.start=start, loc.end= end, num.mark, seg.mean=CN)
write_tsv(IGV_cohort_res, "../ASCAT/IGV_files/IGV_segmsRC_pq_TP53.seg")

