library(ggplot2)
library(plyr)
library(httr)
# BiocManager::install("drawProteins")
library(drawProteins)
library(ggrepel)
library(RODBC)
library(tidyverse)

db <- odbcConnectAccess2007("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/TP53_DB_v1.accdb")

TABS <- sqlTables(db)$TABLE_NAME
TABS

# IMPORT SAMPLES INFO
samples <- sqlFetch(db,"ELENCO_CAMPIONI_TP53_121119")
sampinfo <- samples %>% select(SEQ_ID, DNA, CXFASE)


# IMPORT MUTS
mut <- sqlFetch(db, "All_GOLD_mut_TP53")

# correct annotation on func
mut$ExonicFuncrefGene <- mut$ExonicFuncrefGene %>% str_replace("_"," ")
mut$ExonicFuncrefGene %>% unique

#_____ select transcript and define protein AA change_____

# NM_001126112 canonical isoform

# method 1> specific transcript
mut$AAchange <- mut$AAChangerefGene %>% as.character() %>% str_extract("NM_001126112[0-9A-Za-z:_\\.]+") %>% str_extract("p.*") %>% str_remove("p\\.")

# method 2> first transcript
# mut$AAchange <- sapply(strsplit(mut$AAChangerefGene %>% as.character(), ":"), "[", 5) %>% str_remove("p\\.") %>% str_remove(",TP53")
# mut$transcript <- sapply(strsplit(mut$AAChangerefGene %>% as.character(), ":"), "[", 2)
# mut$transcript %>% table


#______ exploration analysis on transcripts ______
ISOFORMS <- sapply(mut$AAChangerefGene, str_extract_all, "NM_[0-9]+")

all_iso <- ISOFORMS %>% unlist() %>% unique

matIso <- matrix(data = NA, nrow = length(ISOFORMS), ncol = length(all_iso), dimnames = list(mut$ID, all_iso))
matIsoAA <- matrix(data = NA, nrow = length(ISOFORMS), ncol = length(all_iso), dimnames = list(mut$ID, all_iso))

i=1
for(i in 1:length(all_iso)){
  iso=all_iso[i]
  
  j=10
  for (j in 1:length(ISOFORMS)){
    res <- iso %in% ISOFORMS[[j]] %>% as.numeric()
    matIso[j,i] <- res
    
    resAA <- str_extract(mut$AAChangerefGene[j], paste0(iso,"[0-9A-Za-z:_\\.]+")) %>% str_extract("p.*") %>% str_remove("p\\.")
    matIsoAA[j,i] <- resAA
  }
}

allsum <- apply(matIso, 2, sum)

names(allsum) <- all_iso
allsum

matIso %>% pheatmap::pheatmap()

#_______ format for MutPlot _____

# set AA change for noncoding mutations (splicing and UTR)
mut$FuncrefGene[is.na(mut$AAchange)] %>% as.character()
mut$GeneDetailrefGene[is.na(mut$AAchange)] %>% str_extract("NM_001126112[0-9A-Za-z:_\\.]+") # %>% parse_number()

# AA.change
mut$AAchange[is.na(mut$AAchange)] <- mut$GeneDetailrefGene[is.na(mut$AAchange)] %>% 
  str_match("(NM_001126112)(:exon[0-9]:c\\.)([0-9A-Za-z:_\\.]+)") %>% 
  .[,4] %>% 
  parse_number() %>% 
  `/`(3) %>% 
  round

# UTR AAchange
mut$AAchange[mut$FuncrefGene=="UTR5"] <- 0

# func
mut$FuncrefGene[is.na(mut$ExonicFuncrefGene)] %>% as.character()
mut$ExonicFuncrefGene[is.na(mut$ExonicFuncrefGene)] <- mut$FuncrefGene[is.na(mut$ExonicFuncrefGene)] %>% as.character()



template2 <- data.frame( Hugo_Symbol= mut$GenerefGene, 
                         Sample_ID=mut$DNA, 
                         Protein_Change=mut$AAchange, 
                         Mutation_Type=mut$ExonicFuncrefGene )


template3 <- data.frame( Hugo_Symbol= mut$GenerefGene, 
                         Sample_ID=mut$DNA, 
                         Protein_Change=mut$AAchange, 
                         Mutation_Type=mut$ExonicFuncrefGene,
                         callset=mut$dataset)


# templateInfo <- left_join( template2, sampinfo, by=c("Sample_ID"="DNA") )

#============================ LOLLIPOP PLOT ===================================

# BiocManager::install("drawProteins")

library(ggplot2)
library(plyr)
library(httr)
library(drawProteins)
library(ggrepel)

setwd("D:/Cartelle_personali/Andrea/GIT_Shirke019/Packages/Mutplot")

proteins_acc<-read.table("UniProt.txt",header=T)

template<-read.table("upload.data.txt",header=T)

template<-template[,1:4]

      var<-template2
      
      vardup <- var
      vardup$dup <- as.numeric(duplicated(var))
      
      var<-unique(var[!duplicated(var),])# remove duplicates
      
      var[order(var$Hugo_Symbol,gsub("([A-Z]+)([0-9]+)","\\2",var$Protein_Change)),]#order
      
      var.freq<-plyr::count(var,c("Hugo_Symbol","Protein_Change","Mutation_Type"))
      
      var.aanum<-unlist(
        lapply(
          regmatches(
            var.freq$Protein_Change,gregexpr(pattern="*(\\d+)",var.freq$Protein_Change)
            ),function(x) x[[1]]))
      
      var.plot.data<-cbind(var.freq,as.numeric(as.character(var.aanum)))#hard way to remove factor
      
      colnames(var.plot.data)[ncol(var.plot.data)]<-"var.aanum"

    GENE="TP53"
    
    AAfreq=1
    
    var.plot <- var.plot.data
    var.plot<-var.plot[var.plot$Hugo_Symbol==GENE,]
    
    var.plot<-var.plot[order(var.plot$var.aanum),]
    var.plot$Protein_Change<-as.character(var.plot$Protein_Change)
    
    p<-ggplot(var.plot,aes(var.aanum,freq,color=Mutation_Type,label=Protein_Change))+
      geom_segment(aes(x=var.aanum,y=0,xend=var.aanum,yend=freq),
                   color=ifelse(var.plot$freq>= AAfreq,"grey70","grey50"),
                   size=ifelse(var.plot$freq>= AAfreq,0.8,0.4))+
      geom_point(size=5)+
      geom_text_repel(nudge_y=0.2,color=ifelse(var.plot$freq>= AAfreq,"black","NA"),
                      size=ifelse(var.plot$freq>= AAfreq,3,2),fontface="bold")
    
    p
    
    
    p_acc2<-proteins_acc[proteins_acc$Hugo_Symbol==GENE,2]
    
    proteins_acc_url<-gsub(" ","%2C",p_acc2)
    
    baseurl<-"https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession="
    
    url<-paste0(baseurl,proteins_acc_url)
    
    prots_feat<-GET(url,accept_json())
    
    prots_feat_red<-httr::content(prots_feat)
    
    features_total_plot<-NULL
    
    
    for(i in 1:length(prots_feat_red)){ 
      print(i)
      features_temp<-drawProteins::extract_feat_acc(prots_feat_red[[i]])#the extract_feat_acc() function takes features into a data.frame
      features_temp$order<-i # this order is needed for plotting later
      features_total_plot<-rbind(features_total_plot,features_temp)
    }
    
    plot_start<-0#starts 0
    plot_end<-max(features_total_plot$end)
    
    if("CHAIN"%in%(unique(features_total_plot$type))) {
      p <- p +
        geom_rect(data=features_total_plot[features_total_plot$type=="CHAIN",],
                  mapping=aes(xmin=begin,xmax=end,ymin=-0.35,ymax=-0.1),
                  colour="grey",fill="grey",
                  inherit.aes=F)
    }
    
    if("DNA_BIND"%in%(unique(features_total_plot$type))) {
      features_total_plot[features_total_plot$type=="DNA_BIND","description"]<-features_total_plot[features_total_plot$type=="DNA_BIND","type"]
      
      p <- p +
        geom_rect(data=features_total_plot[features_total_plot$type=="DNA_BIND",],
                  mapping=aes(xmin=begin,xmax=end,ymin=-0.45,ymax=0,fill=description),
                  inherit.aes=F)
    }
    
    if("DOMAIN"%in%(unique(features_total_plot$type))) {
      p <- p +
        geom_rect(data=features_total_plot[features_total_plot$type=="DOMAIN",],
                  mapping=aes(xmin=begin,xmax=end,ymin=-0.45,ymax=0,fill=description),
                  inherit.aes=F)
    }
    
    if("MOTIF"%in%(unique(features_total_plot$type))) {
      p <- p +
        geom_rect(data=features_total_plot[features_total_plot$type=="MOTIF",],
                  mapping=aes(xmin=begin,xmax=end,ymin=-0.45,ymax=0,fill=description),
                  inherit.aes=F)
    }
    
    p<- p+
      scale_x_continuous(breaks=round(seq(plot_start,plot_end,by=50)),name="Amino acid number") +
      theme(plot.title=element_text(face="bold",size=(15),hjust=0),axis.text.x=element_text(angle=90,hjust=1)) +
      labs(y="Mutation frequency in our dataset",title=paste("Mutation frequency in ",GENE," - complete dataset",sep=""),
           subtitle="reference isoform: NM_001126112\nprotein structure source:Uniprot") +
      scale_y_continuous(breaks=seq(0,max(var.plot$freq),1),name="frequency")
    
    p+theme_classic()+theme(legend.position="bottom")
    
    
    
    
    
    
    p_acc<-proteins_acc[proteins_acc$Hugo_Symbol==GENE,2]
    proteins_acc_url<-gsub(" ","%2C",p_acc)
    baseurl<-"https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession="
    url<-paste0(baseurl,proteins_acc_url)
    prots_feat<-GET(url,accept_json())
    prots_feat_red<-httr::content(prots_feat)
    features_total_plot<-NULL
    for(i in 1:length(prots_feat_red)){ 
      features_temp<-drawProteins::extract_feat_acc(prots_feat_red[[i]])#the extract_feat_acc() function takes features into a data.frame
      features_temp$order<-i # this order is needed for plotting later
      features_total_plot<-rbind(features_total_plot,features_temp)
    }
    
    features_total_plot<-subset(features_total_plot,type%in% c("CHAIN","DNA_BIND","DOMAIN","MOTIF"))
    
  