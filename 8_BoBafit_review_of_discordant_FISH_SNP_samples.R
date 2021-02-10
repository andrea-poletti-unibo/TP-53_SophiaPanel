library(tidyverse)
library(data.table)


#___ Import the file with the samples names _____
import <- readxl::read_excel("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT BOBaFIT - Documenti/run_090221/sample.xlsx")

SNPS <- import$sample

rawcopydir <- "D:/RAW_DATA/RAWCOPY/alpha=10^-7/" # set the rawcopy folders dir

#___ check if all samples selected are present ____

dirs <- list.dirs(rawcopydir, full.names = F, recursive = F)

SNPS %in% dirs %>% table # check
SNPS[!c(SNPS %in% dirs)]


#____ loop over the Rawcopy dirs to extract the segments ____
results <- data.frame()


for(i in SNPS){
  
  print(i)
  
  setwd(paste0(rawcopydir,"/",i))
  
  segms <- data.table::fread("segments2.txt") 
  names(segms)[5] <- "logR"
  names(segms)[7] <- "Allelic_imbalance"
  
  namedsegms <- cbind(sample=i, segms)
  
  results <- rbind(results,namedsegms)
  
}

#____ save the results ____
readr::write_tsv(results, "C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT BOBaFIT - Documenti/run_090221/new_segments.txt")





#################### PART 2 ##########################


library(data.table)
library(tidyverse)
library(GenomicRanges)


segms <- fread("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT BOBaFIT - Documenti/run_090221/new_segments.txt") # import the ASCAT segments

segms <- segms %>% filter( !grepl("X|Y|M", segms$Chromosome)  )

segms$arm <- segms$Start.band %>% str_extract("[0-9]+[pq]")

write_tsv(segms, "C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT BOBaFIT - Documenti/run_090221/pq_profiles_cohort.txt")


############### PART 3 ###################


library(BOBaFIT)


chrlist <- c("10q","2q","10p","12q","2p","16p","20q","17q","12p", 
             "1p" ,"8q","4q","4p","20p","17p","6p","22q","18p","18q","6q")

seg<- fread("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT BOBaFIT - Documenti/run_090221/pq_profiles_cohort.txt")

seg$Chromosome <- seg$Chromosome %>% str_remove("chr") %>% as.numeric()

seg$width <- seg$End - seg$Start

seg$CN <- 2^((seg$logR/0.55) +1)

seg$CN %>% density %>% plot(xlim=c(0,4))


seg$ID <- seg$sample
seg<-seg %>% select(ID,chr= Chromosome,start=Start,end=End, arm,CN, width)
seg_formatted <- as.data.frame(seg)

path <- "C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT BOBaFIT - Documenti/run_090221/plots/"


results_seg<- DRrefit(segments_chort = seg_formatted, 
                      chrlist =chrlist, 
                      maxCN = 6, 
                      clust_method = "ward.D2", 
                      plot_output = T, 
                      plot_path = path )


debugonce(DRrefit)



report <- results_seg$report

write_tsv(report, "C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT BOBaFIT - Documenti/run_090221/report_sample.txt")
