library(data.table)
library(tidyverse)

segm_R_TP53<-data.table::fread("D:/analisi_in_corso/TP53/CN_calls/rawcopy_segm_TP53_cohort.txt")

segm <- segm_R_TP53 %>% select( -Refseq.genes) 
segm$Chromosome <- segm$Chromosome %>%  str_replace("chr","")

segm$CN <- 2^((segm$logR/0.55)+1)
segm$CN[segm$CN>6] <- 6

plot(density(segm$CN, na.rm=T), xlim=c(0,5))



genesLocations <- data.table::fread("D:/analisi_in_corso/TP53/CN_calls/25_genes_loci.txt")

library(GenomicRanges)

p <- makeGRangesFromDataFrame(segm, keep.extra.columns = T)

p




i=2
results <- data.frame()
for (i in 1:25){
  
  gene <- genesLocations$gene[i]
  chr <- as.numeric(genesLocations$chromosome[i])
  start <- as.numeric(genesLocations$start[i])
  end <- as.numeric(genesLocations$end[i])
  
  print(paste(gene,chr, start,end))
  
  Grange <- GRanges(chr, IRanges(start,end))
  
  overlap <- pintersect( p, Grange)
  res.df <- overlap %>% as.data.frame() %>% filter(hit==TRUE)
  
  res.df <- res.df %>% arrange(sample, desc(abs(CN-2)) )%>% distinct(sample , .keep_all=T) # this will keep only the most exteme CN overlapping segm per sample
  
  #_ _ _ _ _ _
  #
  # hits %>% as.data.frame()
  # subjectHits(hits)
  # 
  # # overlaps <- pintersect(p[queryHits(hits)], Grange[subjectHits(hits)])
  # 
  # overlaps <- pintersect(p[queryHits(hits)], Grange)
  # overlaps
  # 
  # overlaps %>% as.data.frame() %>% View()
  # res %>% as.data.frame() %>% View()
  # 
  # Grange[subjectHits(hits)]
  # sub
  # percentOverlap <- width(overlaps) / width(Grange[subjectHits(hits)])
  # hits <- hits[percentOverlap > 0.5]
  #
  #_ _ _ _ _ _ 
  
  
  
  CNvalues<- res.df %>% dplyr::select(sample,CN)
  CNvalues$gene <- gene
  
  results <- rbind(results,CNvalues)
}  


results.table <- spread(results,"gene", "CN")


write_tsv(results,"D:/analisi_in_corso/TP53/CN_calls/calls/CN_Sophia_25genes_melted.tsv")
write_tsv(results.table,"D:/analisi_in_corso/TP53/CN_calls/calls/CN_Sophia_25genes_table.tsv")



######################### PURITY CORRECTED CN CALLS ###############################################

purity <- fread("D:/analisi_in_corso/TP53/SNP_TP53_purity_reviewed.txt")

res.table.purity <- left_join(results.table, purity, by=c("sample"="SNP_file"))

genes <- genesLocations$gene %>% sort

i=genes[1]

for( i in genes) {
  print(i)
  res.table.purity[,paste0(i,"_adj")] = ((res.table.purity[,i] -2) / (res.table.purity$purity/100)) + 2
  
}

res.table.purity %>% select(BRAF, purity, BRAF_adj) %>% View

write_tsv(res.table.purity,"D:/analisi_in_corso/TP53/CN_calls/calls/CN_Sophia_25genes_table_PURITY_ADJ.tsv")

