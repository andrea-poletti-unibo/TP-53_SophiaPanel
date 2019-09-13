library(data.table)
library(tidyverse)


segms <- fread("../ASCAT/All_segmsRC_pq.txt")

samples <- unique(segms$sample)

chrarmdata <- fread("D:/Cartelle_personali/Andrea/GIT_Shirke019/Utility/chrom_arm_data_v2.txt") # import the chromosome arms data (previously generated file)
chrarmdata$length <- chrarmdata$end - chrarmdata$start

arms <- chrarmdata$chrarm


THRESH <- 0.4
broad_event <- 0.25

#__________ ARM call LOOP ___________

CallAmpMatrix <- matrix(NA, nrow = length(samples), ncol = length(arms), dimnames = list(samples, paste0("amp_",arms)))
CallDelMatrix <- matrix(NA, nrow = length(samples), ncol = length(arms), dimnames = list(samples, paste0("del_",arms)))

ValueMatrix <- matrix(NA, nrow = length(samples), ncol = length(arms), dimnames = list(samples, arms))

i = 1
for(i in 1:length(samples)){
  
  samp.seg <- segms %>% filter(sample == samples[i])
  
  j=1
  for (j in 1:length(arms)){
    
    arm_len <- chrarmdata %>% filter(chrarm==arms[j]) %>% .$length
    
    #_______CALLS_________
    # dels
    arm.del.thresh <- samp.seg %>% filter(arm == arms[j], CN < 2-THRESH)
    del.thresh.len <- arm.del.thresh$width %>% sum
    del.call <- if_else(del.thresh.len > (arm_len * broad_event), 1, 0)
    
    CallDelMatrix[i,j] <- del.call # save
    
    # amps
    arm.amp.thresh <- samp.seg %>% filter(arm == arms[j], CN > 2+THRESH)
    amp.thresh.len <- arm.amp.thresh$width %>% sum
    amp.call <- if_else(amp.thresh.len > (arm_len * broad_event), 1, 0)
    
    CallAmpMatrix[i,j] <- amp.call # save
    
    #_______VALUE_________
    
    arm.seg <- samp.seg %>% filter(arm == arms[j])
    arm.value <- matrixStats::weightedMedian(arm.seg$CN, arm.seg$width)
    
    ValueMatrix[i,j] <- arm.value # save
  }
}

#========== Export results ===========

complete.results <- cbind(sample=rownames(CallDelMatrix),CallDelMatrix, CallAmpMatrix) %>% as.data.frame()
write_tsv(complete.results, "../ASCAT/ARM_CALLS/arm_calls_TP53cohort.txt")

value.results <- cbind(sample=rownames(CallDelMatrix),ValueMatrix) %>% as.data.frame()
write_tsv(value.results, "../ASCAT/ARM_CALLS/arm_values_TP53cohort.txt")

  
vec <- cbind(value.results, complete.results)
write_tsv(vec, "../ASCAT/ARM_CALLS/arm_calls+values_TP53cohort.txt")
