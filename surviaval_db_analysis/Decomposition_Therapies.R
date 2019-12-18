library(RODBC)
library(readxl)
library(tidyverse)



df<-read_xlsx("C:/Users/emat/Documents/ELENCO_PAZIENTI_TP53_121119.xlsx")

df1<- df %>% select(UPN, FIRST_LINE_teraphy_AS_Treated)
 table(df1$FIRST_LINE_teraphy_AS_Treated)
 unique(df1$FIRST_LINE_teraphy_AS_Treated) %>% sort()

 df1$FIRST_LINE_FORMATTED <-gsub("TX_AUTO", "TX", df1$FIRST_LINE_teraphy_AS_Treated) 
 
 df1$FIRST_LINE_FORMATTED <-gsub("mant|cons", "", df1$FIRST_LINE_FORMATTED)
 
 df1$FIRST_LINE_FORMATTED <-gsub("HDM-1", "TX", df1$FIRST_LINE_FORMATTED) 
 
 df1$FIRST_LINE_FORMATTED <-gsub("HDM-2", "2TX", df1$FIRST_LINE_FORMATTED) 
 

 unique(df1$FIRST_LINE_FORMATTED) %>% sort()
df1$TX<- ifelse(grepl("TX",df1$FIRST_LINE_FORMATTED),1,0)

df1$TX_2<- ifelse(grepl("2TX",df1$FIRST_LINE_FORMATTED),1,0)

df1$BORTEZOMIB<- ifelse(grepl("VCD|BORTEZOMIB|VMP|VD", df1$FIRST_LINE_FORMATTED),1,0)

df1$LENALIDOMIDE<- ifelse(grepl("RD|LENA|REV|KRD|KRd|KR", df1$FIRST_LINE_FORMATTED),1,0)

df1$TALIDOMIDE<- ifelse(grepl("TD|TALIDOMIDE|VTD", df1$FIRST_LINE_FORMATTED),1,0)

df1$CICLOFOSFAMIDE<- ifelse(grepl("KCD|VCD|CRD", df1$FIRST_LINE_FORMATTED),1,0)

df1$DESAMETASONE<- ifelse(grepl("VTD|DEX|dex|KRD|KCD|RD|TD|VD", df1$FIRST_LINE_FORMATTED),1,0)

df1$INTERFERONE<- ifelse(grepl("IFN", df1$FIRST_LINE_FORMATTED),1,0)

df1$CORTISONE<- ifelse(grepl("CORTISONE", df1$FIRST_LINE_FORMATTED),1,0)

df1$PREDNISONE<- ifelse(grepl("VMP|MP|MPV", df1$FIRST_LINE_FORMATTED),1,0)

df1$MELPHALAN<- ifelse(grepl("VMP|MPV|MEL200|MP", df1$FIRST_LINE_FORMATTED),1,0)

df1$VINCRISTINA<- ifelse(grepl("VAD", df1$FIRST_LINE_FORMATTED),1,0)

df1$DEXORUBICINA<- ifelse(grepl("VAD", df1$FIRST_LINE_FORMATTED),1,0)

df1$PEMBROLIZUMAB<- ifelse(grepl("PEMBRO", df1$FIRST_LINE_FORMATTED),1,0)

df1$CARFIZOMIB<-ifelse(grepl("KRD|KRd|KR|KCD", df1$FIRST_LINE_FORMATTED),1,0)

df1$IMID_CLASS<- ifelse(df1$LENALIDOMIDE==1|df1$TALIDOMIDE==1,1,0)
df1$PI_CLASS<-ifelse(df1$BORTEZOMIB==1|df1$CARFIZOMIB==1,1,0)
df1$CORTICOSTEROID_CLASS<-ifelse(df1$DESAMETASONE==1|df1$CORTISONE==1|df1$PREDNISONE==1,1,0)


write_tsv(df1, "C:/Users/emat/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/FL_Therapies_table.txt")




