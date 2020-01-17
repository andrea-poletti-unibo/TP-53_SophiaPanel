#############################################################################
########################### TP53 PFS2 analysis ##############################
#############################################################################

library(RODBC)
library(tidyverse)

db <- odbcConnectAccess2007("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/TP53_DB_v1.accdb")

TABS <- sqlTables(db)$TABLE_NAME
TABS

# quesry per tutti i pazienti del DB che hanno i dati molecolari di relapse (65)
query_relapse <- sqlFetch(db, "Copia di Query_Survival_Analysis")

# filtro per i pazienti che sono progrediti (2 eccezioni)
prog_true <- query_relapse %>% filter(Prog_I==1)


# completamento della colonna prog_II: inserimento di 0 nei valori mancanti (pazienti che non sono progrediti ma sono morti o ancora vivi, mancava il dato)
prog_true$Prog_II[is.na(prog_true$Prog_II)] <- 0


# calcolo nuova variabile = PFS_2_event
prog_true$PFS_2_event <- ifelse( prog_true$Prog_II == 1 | prog_true$OS_event_death == 1, 1, 0)

prog_true$PFS_2_event %>% table

# convert the dates in correct Date format
prog_true$data_II_relapse <- as.Date(prog_true$data_II_relapse, origin="1970-01-01")
prog_true$data_II_relapse %>% class

prog_true$Date_of_death <- as.Date(prog_true$Date_of_death, origin="1970-01-01")
prog_true$Date_of_death %>% class

prog_true$LAST_FULLOW_UP <- as.Date(prog_true$LAST_FULLOW_UP, origin="1970-01-01")
prog_true$LAST_FULLOW_UP %>% class

prog_true$PFS_I_date <- as.Date(prog_true$PFS_I_date, origin="1970-01-01")
prog_true$PFS_I_date %>% class

# calcolo nuova variabile = PFS_2_date

prog_true$PFS_2_date <- dplyr::if_else(prog_true$Prog_II == 1, 
                                       as.Date(prog_true$data_II_relapse, origin="1970-01-01"), 
                                       dplyr::if_else(prog_true$OS_event_death == 1, 
                                                      as.Date(prog_true$Date_of_death, origin="1970-01-01"), 
                                                      as.Date(prog_true$LAST_FULLOW_UP,origin="1970-01-01")
                                       ))

prog_true %>% select(Prog_II, data_II_relapse, OS_event_death, Date_of_death, LAST_FULLOW_UP, PFS_2_event, PFS_2_date) %>% View


# calulate the PFS_2_time

prog_true$PFS_2_time <- prog_true$PFS_2_date - prog_true$PFS_I_date
prog_true %>% select( PFS_I_date, PFS_2_date, PFS_2_time) %>% View

# calulate the PFS_2_months

prog_true$PFS_2_months <- round(prog_true$PFS_2_time / 30.5, 0) 
prog_true %>% select( PFS_I_date, PFS_2_date, PFS_2_time, PFS_2_months) %>% View


################################## survival analysys ############################################################

library(data.table)
library(tidyverse)
library(survival)
library(survminer)
library(clusteval)

prog_true$flag <- ifelse(prog_true$PROTOCOLLO != "BO2005" & prog_true$PROTOCOLLO != "EMN02" & 
                           prog_true$TTP_I_months <= 24 & prog_true$Prog_I==1 &
                           prog_true$TP53_adj_D >= 1.9 & prog_true$MUT_p53_D_SUB_CLON ==0 , 
                         1, 0)

prog_true$flag %>% table

DF2 <- prog_true %>% filter(flag==0)


# ISS distribution Check
DF2$ISS %>% table 


# del 10%
DF2$call10_del_TP53_D <- ifelse(DF2$TP53_adj_D <= 1.9 , 1, 0)
DF2$call10_del_TP53_R <- ifelse(DF2$TP53_adj_R <= 1.9 , 1, 0)

DF2$call10_del_TP53_D %>% table
DF2$call10_del_TP53_R %>% table

DF2 %>% select(call10_del_TP53_D, call10_del_TP53_R, TP53_adj_D, TP53_adj_R) %>% View

table(DF2$call10_del_TP53_D, DF2$call10_del_TP53_R, useNA = c("always"))

cor.test(DF2$call10_del_TP53_D, DF2$call10_del_TP53_R)
cluster_similarity(DF2$call10_del_TP53_D, DF2$call10_del_TP53_R, similarity="jaccard", method="independence")


# del 50%
DF2$call50_del_TP53_D <- ifelse(DF2$TP53_adj_D <= 1.5 , 1, 0)
DF2$call50_del_TP53_R <- ifelse(DF2$TP53_adj_R <= 1.5 , 1, 0)

DF2$call50_del_TP53_D %>% table
DF2$call50_del_TP53_R %>% table

DF2 %>% select(call10_del_TP53_D, call10_del_TP53_R, TP53_adj_D, TP53_adj_R) %>% View

table(DF2$call50_del_TP53_D, DF2$call50_del_TP53_R, useNA = c("always"))

cor.test(DF2$call50_del_TP53_D, DF2$call50_del_TP53_R)
cluster_similarity(DF2$call50_del_TP53_D, DF2$call50_del_TP53_R, similarity="jaccard", method="independence")


#______________________
# del 10% survival PFS2
#______________________

PFS2 <- Surv( time = DF2$PFS_2_months, event = DF2$PFS_2_event)

# DEL at RELAPSE
ggsurvplot(survfit(PFS2 ~ DF2$call10_del_TP53_R, data = DF2), pval = T, risk.table = T, xlab = "PFS")

coxph(PFS2 ~ DF2$call10_del_TP53_R, data = DF2) %>% summary
coxph(PFS2 ~ DF2$call10_del_TP53_R + strata(DF2$ISS) , data = DF2) %>% summary

# DEL at DIAGNOSIS
ggsurvplot(survfit(PFS2 ~ DF2$call10_del_TP53_D, data = DF2), pval = T, risk.table = T, xlab = "PFS")

coxph(PFS2 ~ DF2$call10_del_TP53_D, data = DF2) %>% summary
coxph(PFS2 ~ DF2$call10_del_TP53_D + strata(DF2$ISS) , data = DF2) %>% summary

#______________________
# mut survival PFS2
#______________________

cor.test(DF2$MUT_p53_D_SUB_CLON, DF2$MUT_P53_R_SUB_CLON)
cluster_similarity(DF2$MUT_p53_D_SUB_CLON, DF2$MUT_P53_R_SUB_CLON, similarity="jaccard", method="independence")

# MUT at RELAPSE
ggsurvplot(survfit(PFS2 ~ DF2$MUT_P53_R_SUB_CLON, data = DF2), pval = T, risk.table = T, xlab = "PFS")

coxph(PFS2 ~ DF2$MUT_P53_R_SUB_CLON, data = DF2) %>% summary
coxph(PFS2 ~ DF2$MUT_P53_R_SUB_CLON + strata(DF2$ISS) , data = DF2) %>% summary

# MUT at DIAGNOSIS
ggsurvplot(survfit(PFS2 ~ DF2$MUT_p53_D_SUB_CLON, data = DF2), pval = T, risk.table = T, xlab = "PFS")

coxph(PFS2 ~ DF2$MUT_p53_D_SUB_CLON, data = DF2) %>% summary
coxph(PFS2 ~ DF2$MUT_p53_D_SUB_CLON + strata(DF2$ISS) , data = DF2) %>% summary

#______________________
# mut-del survival PFS2
#______________________

DF2$mut_del_R <- ifelse(DF2$call10_del_TP53_R==1 & DF2$MUT_P53_R_SUB_CLON ==1, 1, 0)
DF2$mut_del_R %>% table

DF2$mut_del_D <- ifelse(DF2$call10_del_TP53_D==1 & DF2$MUT_p53_D_SUB_CLON ==1, 1, 0)
DF2$mut_del_D %>% table

table(DF2$mut_del_D, DF2$mut_del_R)

# MUT-DEL at RELAPSE
ggsurvplot(survfit(PFS2 ~ DF2$mut_del_R, data = DF2), pval = T, risk.table = T, xlab = "PFS")

coxph(PFS2 ~ DF2$mut_del_R, data = DF2) %>% summary
coxph(PFS2 ~ DF2$mut_del_R + strata(DF2$ISS) , data = DF2) %>% summary

# MUT-DEL at DIAGNOSIS
ggsurvplot(survfit(PFS2 ~ DF2$mut_del_D, data = DF2), pval = T, risk.table = T, xlab = "PFS")

coxph(PFS2 ~ DF2$MUT_p53_D_SUB_CLON, data = DF2) %>% summary
coxph(PFS2 ~ DF2$MUT_p53_D_SUB_CLON + strata(DF2$ISS) , data = DF2) %>% summary
