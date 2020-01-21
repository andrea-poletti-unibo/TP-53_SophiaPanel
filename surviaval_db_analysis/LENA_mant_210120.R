library(RODBC)
library(tidyverse)
library(data.table)
library(survival)
library(survminer)

################## data import and PFS 2 creation ################
db <- odbcConnectAccess2007("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/TP53_DB_v1.accdb")

query_relapse <- sqlFetch(db, "Copia di Query_Survival_Analysis")
prog_true <- query_relapse %>% filter(Prog_I==1)
prog_true$Prog_II[is.na(prog_true$Prog_II)] <- 0

# calcolo nuova variabile = PFS_2_event
prog_true$PFS_2_event <- ifelse( prog_true$Prog_II == 1 | prog_true$OS_event_death == 1, 1, 0)

# convert the dates in correct Date format
prog_true$data_II_relapse <- as.Date(prog_true$data_II_relapse, origin="1970-01-01")
prog_true$Date_of_death <- as.Date(prog_true$Date_of_death, origin="1970-01-01")
prog_true$LAST_FULLOW_UP <- as.Date(prog_true$LAST_FULLOW_UP, origin="1970-01-01")
prog_true$PFS_I_date <- as.Date(prog_true$PFS_I_date, origin="1970-01-01")

# calcolo nuova variabile = PFS_2_date
prog_true$PFS_2_date <- dplyr::if_else(prog_true$Prog_II == 1, 
                                       as.Date(prog_true$data_II_relapse, origin="1970-01-01"), 
                                       dplyr::if_else(prog_true$OS_event_death == 1, 
                                                      as.Date(prog_true$Date_of_death, origin="1970-01-01"), 
                                                      as.Date(prog_true$LAST_FULLOW_UP,origin="1970-01-01")
                                       ))

# calulate the PFS_2_time
prog_true$PFS_2_time <- prog_true$PFS_2_date - prog_true$PFS_I_date
# calulate the PFS_2_months
prog_true$PFS_2_months <- round(prog_true$PFS_2_time / 30.5, 0) 


# FLAG filter
prog_true$flag <- ifelse(prog_true$PROTOCOLLO != "BO2005" & prog_true$PROTOCOLLO != "EMN02" & 
                           prog_true$TTP_I_months <= 24 & prog_true$Prog_I==1 &
                           prog_true$TP53_adj_D >= 1.9 & prog_true$MUT_p53_D_SUB_CLON ==0 , 
                         1, 0)

DF2 <- prog_true %>% filter(flag==0)

# calulate the PFS2 object
PFS2 <- Surv( time = DF2$PFS_2_months, event = DF2$PFS_2_event)

# call del 10% creation
DF2$call10_del_TP53_D <- ifelse(DF2$TP53_adj_D <= 1.9 , 1, 0)
DF2$call10_del_TP53_R <- ifelse(DF2$TP53_adj_R <= 1.9 , 1, 0)



####################### distribuzione della TTP1 nei pts con progessione ########################

dels <- DF2 %>% filter(DF2$call10_del_TP53_R==1)
nodels <- DF2 %>% filter(DF2$call10_del_TP53_R==0)


dels$TTP_I_months %>% sort
nodels$TTP_I_months %>% sort

dels$TTP_I_months %>% summary
nodels$TTP_I_months %>% summary


t.test(dels$TTP_I_months, nodels$TTP_I_months) 
# i due gruppi non hanno un distribuzione di tempo di progressione 1 significativamente diversa! No necessità di correggere per differenze di TTP1


####################### MANTENIMENTO con LENA ########################

table(DF2$MAINTENANCE_yes_1_no_0, DF2$call10_del_TP53_R)

ggsurvplot(survfit(PFS2 ~ DF2$MAINTENANCE_yes_1_no_0, data = DF2), pval = T, risk.table = T, xlab = "PFS")

DF2$Maintenance_Therapy %>% table

DF2$LENA_MANT <- ifelse(DF2$Maintenance_Therapy == "LENA", 1, 0)

DF2$LENA_MANT %>% table

# il mantenimento con lena è significativo? Qui NO.... ma da altri lavori si sà che lo è!
ggsurvplot(survfit(PFS2 ~ DF2$LENA_MANT, data = DF2), pval = T, risk.table = T, xlab = "PFS")
coxph(PFS2 ~ DF2$LENA_MANT, data = DF2) %>% summary

gmodels::CrossTable(DF2$LENA_MANT, DF2$call10_del_TP53_R, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F, fisher = T)


# i pts con mantenimento con LENA come vanno in TTP2 in base al fatto che siano deleti o meno su TP53?
LENAmant_DF2 <- DF2 %>% filter(DF2$LENA_MANT==1)

PFS2_lena <- Surv( time = LENAmant_DF2$PFS_2_months, event = LENAmant_DF2$PFS_2_event)

ggsurvplot(survfit(PFS2_lena ~ LENAmant_DF2$call10_del_TP53_R, data = LENAmant_DF2), pval = T, risk.table = T, xlab = "PFS")
coxph(PFS2_lena ~ LENAmant_DF2$call10_del_TP53_R, data = LENAmant_DF2) %>% summary
