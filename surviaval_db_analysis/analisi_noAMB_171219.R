library(tidyverse)
library(RODBC)

db <- odbcConnectAccess2007("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/TP53_DB_v1-AndreaSB.accdb")

df <- sqlFetch(db, "Query_Survival_Analysis")

df$PROTOCOLLO %>% table

amb <- df %>% filter(PROTOCOLLO != "BO2005" & PROTOCOLLO != "EMN02")

amb$PROTOCOLLO %>% table

amb$FIRST_LINE_FORMATTED %>% table

amb12 <- amb %>% filter(TTP_I_months <= 12 & Prog_I==1 )

amb12$FIRST_LINE_FORMATTED %>% as.character() %>% table
amb12$TERAPIA_INDUZIONE_I_LINEA %>% as.character() %>% table




amb24 <- amb %>% filter(TTP_I_months <= 24 & Prog_I==1 )

amb24$FIRST_LINE_FORMATTED %>% as.character() %>% table
amb24$TERAPIA_INDUZIONE_I_LINEA %>% as.character() %>% table


table(amb24$TP53_adj < 1.9, amb24$MUT_p53_D_SUB_CLON)


df$flag <- ifelse(df$PROTOCOLLO != "BO2005" & df$PROTOCOLLO != "EMN02" & 
                    df$TTP_I_months <= 24 & df$Prog_I==1 &
                    df$TP53_adj >= 1.9 & df$MUT_p53_D_SUB_CLON ==0 , 
                  1, 0)


table(df$flag)


DF <- df %>% filter(flag != 1)


DF$call10_del_TP53 <- ifelse(DF$TP53_adj <= 1.9, 1,0)
table(DF$call10_del_TP53, DF$Fish_DEL_17)

#============ Survival analysis =====================

library(data.table)
library(tidyverse)
library(survival)
library(survminer)

OS <- Surv( time = DF$OS_MESI, event = DF$OS_event_death)
PFS <- Surv( time = DF$PFS_I_months, event = DF$PFS_I_event)


ggsurvplot(survfit(OS ~ DF$call10_del_TP53, data = DF), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ DF$call10_del_TP53, data = DF), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ DF$call10_del_TP53, data = DF) %>% summary
coxph(PFS ~ DF$call10_del_TP53, data = DF) %>% summary

coxph(OS ~ DF$call10_del_TP53 + DF$ISS, data = DF) %>% summary
coxph(PFS ~ DF$call10_del_TP53 + DF$ISS , data = DF) %>% summary

coxph(OS ~ DF$call10_del_TP53 + strata(DF$ISS), data = DF) %>% summary
coxph(PFS ~ DF$call10_del_TP53 + strata(DF$ISS) , data = DF) %>% summary


gmodels::CrossTable(DF$ISS, DF$call10_del_TP53, fisher = T)
gmodels::CrossTable(DF$ISS==3, DF$call10_del_TP53, fisher = T)
