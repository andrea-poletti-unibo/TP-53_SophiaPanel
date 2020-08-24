library(RODBC)
library(data.table)
library(tidyverse)

db <- odbcConnectAccess2007("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/TP53_DB_v1.accdb")

TABS <- sqlTables(db)$TABLE_NAME
TABS

query_diagnosis <- sqlFetch(db, "Query_Survival_Analysis")


query_diagnosis$flag <- ifelse(query_diagnosis$PROTOCOLLO != "BO2005" & query_diagnosis$PROTOCOLLO != "EMN02" & 
                                 query_diagnosis$TTP_I_months <= 24 & query_diagnosis$Prog_I==1 &
                                 query_diagnosis$TP53_adj >= 1.9 & query_diagnosis$MUT_p53_D_SUB_CLON ==0 , 
                               1, 0)

query_diagnosis$flag %>% table

DF1 <- query_diagnosis %>% filter(flag==0)

# del 10%
DF1$call10_del_TP53_D <- ifelse(DF1$TP53_adj <= 1.9 , 1, 0)

DF1$call10_del_TP53_D %>% table


#______fish snp comparison_____

SNP_call10_del_TP53_D <- DF1$call10_del_TP53_D
Fish_DEL_17p <- DF1$Fish_DEL_17
Fish_DEL_17p[is.na(Fish_DEL_17p)] <- "NA"

gmodels::CrossTable(SNP_call10_del_TP53_D, Fish_DEL_17p, prop.r = F, prop.c = F, prop.chisq = F)

table(DF1$call10_del_TP53_D, DF1$Fish_DEL_17, useNA = "always")

DF1$comparison_F_S <- ifelse(DF1$call10_del_TP53_D==1 & DF1$Fish_DEL_17==1, "pos_concord",
                      ifelse(DF1$call10_del_TP53_D==1 & DF1$Fish_DEL_17==0, "only_SNP",
                      ifelse(DF1$call10_del_TP53_D==0 & DF1$Fish_DEL_17==1, "only_FISH", "neg_concord"))) 

DF1$comparison_F_S %>% table

DF1 %>% filter(comparison_F_S=="only_SNP") %>% select(UPN,CEL_NAME_D,TP53_adj, call10_del_TP53_D,Fish_DEL_17,purity, purity_source) %>% View

DF1 %>% filter(comparison_F_S=="only_SNP") %>% select(UPN,CEL_NAME_D,TP53_adj, call10_del_TP53_D,Fish_DEL_17,purity, purity_source) %>% write.table("clipboard", sep="\t", row.names=FALSE)



