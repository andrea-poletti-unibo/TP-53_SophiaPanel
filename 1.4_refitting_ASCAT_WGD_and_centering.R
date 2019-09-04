library(data.table)
library(tidyverse)

segms <- fread("../ASCAT/All_adj_fitted_and_raw_segments.tsv")

segms_Wmed <- segms %>% group_by(sample) %>% mutate( WeightedMedian= matrixStats::weightedMedian(CNvalueRaw, length))

WM <- segms_Wmed %>% summarise(WM=mean(WeightedMedian))

# plot(WM$WM, col=ifelse(WM$WM>3.5,"red","blue"))
# abline(h=3.5)

#=============== refitting WGD ======================
segms_Wmed$CNvalueRaw_WGDref <- ifelse(segms_Wmed$WeightedMedian > 3.5, segms_Wmed$CNvalueRaw/2, segms_Wmed$CNvalueRaw) 

#=============== centering WGD ====================

# preliminary step for centering: compute a second median on the refitted CN
segms_Wmed2 <- segms_Wmed %>% group_by(sample) %>% mutate( WeightedMedian2= matrixStats::weightedMedian(CNvalueRaw_WGDref, length))

# CENTERING: subtract the 2nd median from the CN, plus 2
segms_Wmed2$CNvalueRaw_WGDref_cen <- (segms_Wmed2$CNvalueRaw_WGDref - segms_Wmed2$WeightedMedian2) +2

# correction: if less then 0 CN set to 0 CN
segms_Wmed2$CNvalueRaw_WGDref_cen <- ifelse(segms_Wmed2$CNvalueRaw_WGDref_cen < 0 , 0, segms_Wmed2$CNvalueRaw_WGDref_cen)


# ================= EXPORT ===================

write_tsv(segms_Wmed2, "../ASCAT/All_segms_refit_center.txt")
