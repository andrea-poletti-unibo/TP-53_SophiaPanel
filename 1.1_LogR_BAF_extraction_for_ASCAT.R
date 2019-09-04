library(tidyverse)
library(data.table)

samples <- readxl::read_xlsx("../ELENCO_CAMPIONI_TP53_definitivo_220719.xlsx") 

files <- samples$CEL_NAME %>% na.omit() %>% unique
files
#====================== define the batches =======================

batch1_6.0 <- files[grep("GenomeWide",files)]
batch2_HD <- files[74:167]
batch3_HD <- files[168:262]

batches <- list(batch1_6.0=batch1_6.0, 
                batch2_HD=batch2_HD, 
                batch3_HD=batch3_HD)

# check
Rawcopy <- list.dirs("D:/RAW_DATA/RAWCOPY/alpha=10^-7/", recursive = F, full.names = F)
c(batch1_6.0, batch2_HD, batch3_HD) %in% Rawcopy %>% table


#___________ Loop for extraction from Rawcopy .Rdata files in batches ______________

project="TP53"

j=1
for(j in 1:length(batches)){
  
  bc <- batches[[j]]
  bc.name <- names(batches[j])
  
  # set base working directory (Rawcopy folders dir)
  basewd<-"D:/RAW_DATA/RAWCOPY/alpha=10^-7/"
  
  # LOOP to extract logR and BAF from the .Rdata rawcopy output of each sample
  result.logR<-list()
  result.BAF<-list()
  
  for (i in bc) { 
    t1<-Sys.time()
    
    print(i)
    setwd(paste0(basewd,"/",i)) # change wd to a new sample
    load("rawcopy.Rdata") #load the sample specific .Rdata file
    result.logR[[i]]<-probes.txt$Value #add the sample logR to the total list
    result.BAF[[i]]<-probes.txt$B.Allele.Frequency # add the sample BAF to the total list
    
    t2<-Sys.time()
    print(t2-t1)
  }
  
  # transforming the list to a data frame
  final.logR<-as.data.frame(result.logR) 
  final.BAF<-as.data.frame(result.BAF)
  
  # add columns with probe info
  complete.logR<-base::cbind(probes.txt[,1:3],final.logR)
  rm(final.logR)
  
  complete.BAF<- base::cbind(probes.txt[,1:3],final.BAF)
  rm(final.BAF)
  
  # removing "chr" in Chromosome column
  complete.logR$Chromosome<- gsub("chr","", complete.logR$Chromosome)
  complete.BAF$Chromosome<- gsub("chr","", complete.BAF$Chromosome)
  
  # renaming columns in ASCAT format
  colnames(complete.logR)[1:3]<-c("probe","chrs","pos")
  colnames(complete.BAF)[1:3]<-c("probe","chrs","pos")
  
  # remove sex chromosomes and Mitochondrial 
  LOGR_noXYM <- filter(complete.logR, chrs != "M",chrs != "X", chrs !="Y")
  BAF_noXYM <- filter(complete.BAF, chrs != "M",chrs != "X", chrs !="Y")
  
  
  # exporting the files
  
  dir.create(paste0("D:/analisi_in_corso/", project, "/ASCAT"))
  
  dir.create(paste0("D:/analisi_in_corso/", project, "/ASCAT/", bc.name))
  
  setwd(paste0("D:/analisi_in_corso/", project, "/ASCAT/", bc.name))
  
  write_tsv( LOGR_noXYM, paste0("ascatReady.noXYM_logR_", project, "_", bc.name,".txt"))
  write_tsv( BAF_noXYM, paste0("ascatReady.noXYM_BAF_", project, "_", bc.name,".txt"))
  
}
