library(data.table)
library(tidyverse)

import <- fread("D:/analisi_in_corso/TP53/CN_calls/rawcopy_segm_TP53_cohort.txt")

samples <- unique(import$sample)

TP53 <- import %>% filter(grepl("TP53, |TP53$",import$Refseq.genes ))

TP53_def <- TP53 %>% filter(!is.na(Allelic_imbalance))
TP53_def <- TP53_def %>% arrange(sample, SNPs) %>% distinct(sample, .keep_all = T)
TP53_fail <- TP53 %>% filter(is.na(Allelic_imbalance))


called_samples <- unique(TP53_def$sample)


samples %in% called_samples %>% table
samples[!c(samples %in% called_samples)]

TP_53_LOH <- TP53_def %>% select(sample:SNPs) %>% mutate(LOH_call= as.numeric(Allelic_imbalance > 0.25 & logR > -0.1 & logR < 0.1) )

write_tsv(TP_53_LOH, "../LOH_calls.txt")


