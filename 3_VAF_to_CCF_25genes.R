#___________________ from Jonas ____________________________
#' Mutation burden to mutation copy number
#' 
#' Function to convert mutation burdens into mutation copy number
#' @param burden A vector containing mutation burdens
#' @param totalCopyNumber A vector with total tumour copynumber
#' @param cellularity Float with the cellularity of this sample
#' @param normalCopyNumber A vector with the total copy number of the normal
#' @return A vector with mutation copy number
#' @author dw9
#' @export

mutationBurdenToMutationCopyNumber = function(burden, totalCopyNumber, cellularity, normalCopyNumber=rep(2, length(burden))) {
  mutCopyNumber = burden/cellularity*(cellularity*totalCopyNumber+normalCopyNumber*(1-cellularity))
  mutCopyNumber[is.nan(mutCopyNumber)] = 0
  return(mutCopyNumber)
}

#' Mutation copy number to mutation burden
#' 
#' Function to convert mutation copy number to mutation burden
#' @param copyNumber A vector containing mutation copy number
#' @param totalCopyNumber A vector with total tumour copynumber
#' @param cellularity Float with the cellularity of this sample
#' @param normalCopyNumber A vector with the total copy number of the normal
#' @return A vector with mutation burdens
#' @author dw9
#' @export

mutationCopyNumberToMutationBurden = function(copyNumber, totalCopyNumber, cellularity, normalCopyNumber=rep(2, length(copyNumber))){
  burden = copyNumber*cellularity/(cellularity*totalCopyNumber+normalCopyNumber*(1-cellularity))
  burden[is.nan(burden) | (burden<0.000001)] = 0.000001
  burden[burden>0.999999] = 0.999999
  return(burden)	
}

###############################################




library(tidyverse)
library(data.table)



database <- fread("D:/analisi_in_corso/TP53/ELENCO_CAMPIONI_TP53_definitivo_220719.txt")

database.info <- database %>% select( CEL_NAME, SEQ_ID, SNPnote= SNP)    #DNA_mut, file, PROJECT_CLONAL_SNP)

purityTab<-fread("D:/analisi_in_corso/TP53/SNP_TP53_purity_reviewed.txt") #import ASCAT samples purity values

CNtable<- fread("D:/analisi_in_corso/TP53/CN_calls/calls/CN_Sophia_25genes_melted.tsv") #import copy number of TP53




#================== Gold_OVER_5 CCF calculation ==================


import_over5 <- readxl::read_excel("D:/analisi_in_corso/TP53/GOLD.Converted_all_categories_over5vaf.xlsx")

import_over5_SNP <- left_join(import_over5, database.info, by=c("SEQ_ID"="SEQ_ID")) 

#CHECK missing snp 
import_over5_SNP %>% filter(is.na(CEL_NAME)) %>% select(conversion, SNPnote) %>% distinct() %>% View


# ________ ADD CN VALUES _________
complete_over5 <- left_join(import_over5_SNP, CNtable, by=c("CEL_NAME"="sample", "Gene.refGene"="gene"))
complete_over5 %>% select(CN,Gene.refGene,CEL_NAME) %>% View

# ______ ADD PURITY VALUES _________
complete_over5 <- left_join(complete_over5, purityTab[,-c(1,5)],  by=c("CEL_NAME"="SNP_file"))
complete_over5$purity <- complete_over5$purity/100

#_______ mutCN (u paramerer) computing _______# MM EQ.3

complete_over5$mutCN <- (complete_over5$VAF / complete_over5$purity) * complete_over5$CN  

complete_over5 %>% select(SEQ_ID, CEL_NAME, Gene.refGene, VAF, CN, purity, mutCN) %>% View
plot(complete_over5$mutCN)

#_______ CCF computing _______# MM EQ.4

complete_over5$CCF <- ifelse(complete_over5$mutCN >= 1, # if mutCN is more then 1 the mut is CLONAL
                             1, # then CCF = 1
                             complete_over5$mutCN) # if mutCN is less then 1 the mut is SUBCLONAL

complete_over5 %>% select(SEQ_ID, CEL_NAME, Gene.refGene, VAF, CN, purity, mutCN, CCF) %>% View

# EXPORT 
write_tsv(complete_over5, "D:/analisi_in_corso/TP53/CCF_results/GOLD.Converted_all_categories_over5_CCF.txt")



#================== Gold_UNDER_5 CCF calculation ==================


import_under5 <- readxl::read_excel("D:/analisi_in_corso/TP53/GOLD.Converted_all_categories_under5vaf.xlsx")

import_under5_SNP <- left_join(import_under5, database.info, by=c("SEQ_ID"="SEQ_ID")) 

#CHECK missing snp 
import_under5_SNP %>% filter(is.na(CEL_NAME)) %>% select(conversion, SNPnote) %>% distinct() %>% View


# ________ ADD CN VALUES _________
complete_under5 <- left_join(import_under5_SNP, CNtable, by=c("CEL_NAME"="sample", "Gene.refGene"="gene"))
complete_under5 %>% select(CN,Gene.refGene,CEL_NAME) %>% View

# ______ ADD PURITY VALUES _________
complete_under5 <- left_join(complete_under5, purityTab[,-c(1,5)],  by=c("CEL_NAME"="SNP_file"))
complete_under5$purity <- complete_under5$purity/100

#_______ mutCN (u paramerer) computing _______# MM EQ.3

complete_under5$mutCN <- (complete_under5$VAF / complete_under5$purity) * complete_under5$CN  

complete_under5 %>% select(SEQ_ID, CEL_NAME, Gene.refGene, VAF, CN, purity, mutCN) %>% View
plot(complete_under5$mutCN)

#_______ CCF computing _______# MM EQ.4

complete_under5$CCF <- ifelse(complete_under5$mutCN >= 1, # if mutCN is more then 1 the mut is CLONAL
                             1, # then CCF = 1
                             complete_under5$mutCN) # if mutCN is less then 1 the mut is SUBCLONAL

complete_under5 %>% select(SEQ_ID, CEL_NAME, Gene.refGene, VAF, CN, purity, mutCN, CCF) %>% View

# EXPORT 
write_tsv(complete_under5, "D:/analisi_in_corso/TP53/CCF_results/GOLD.Converted_all_categories_under5_CCF.txt")





#================== Gold_ROCHE454 CCF calculation ==================


import_roche <- fread("D:/analisi_in_corso/TP53/GOLD_Roche454_Category_A-B.txt")
import_roche$VAF <- import_roche$freq


import_roche_SNP <- left_join(import_roche, database.info, by=c("SEQ_ID"="SEQ_ID")) 

#CHECK missing snp 
import_roche_SNP %>% filter(is.na(CEL_NAME)) %>% select(SEQ_ID, SNPnote) %>% distinct() %>% View


# ________ ADD CN VALUES _________
complete_roche <- left_join(import_roche_SNP, CNtable, by=c("CEL_NAME"="sample", "Gene.refGene"="gene"))
complete_roche %>% select(CN,Gene.refGene,CEL_NAME) %>% View

# ______ ADD PURITY VALUES _________
complete_roche <- left_join(complete_roche, purityTab[,-c(1,5)],  by=c("CEL_NAME"="SNP_file"))
complete_roche$purity <- complete_roche$purity/100

#_______ mutCN (u paramerer) computing _______# MM EQ.3

complete_roche$mutCN <- (complete_roche$VAF / complete_roche$purity) * complete_roche$CN  

complete_roche %>% select(SEQ_ID, CEL_NAME, Gene.refGene, VAF, CN, purity, mutCN) %>% View
plot(complete_roche$mutCN)

#_______ CCF computing _______# MM EQ.4

complete_roche$CCF <- ifelse(complete_roche$mutCN >= 1, # if mutCN is more then 1 the mut is CLONAL
                             1, # then CCF = 1
                             complete_roche$mutCN) # if mutCN is less then 1 the mut is SUBCLONAL

complete_roche %>% select(SEQ_ID, CEL_NAME, Gene.refGene, VAF, CN, purity, mutCN, CCF) %>% View

# EXPORT 
write_tsv(complete_roche, "D:/analisi_in_corso/TP53/CCF_results/GOLD_roche454_categories_A-B_CCF.txt")


