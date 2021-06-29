library(RODBC)
library(data.table)
library(tidyverse)
library(survival)
library(survminer)
library(broom)

# full file path to Access DB
file_path <- "C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/TP53_DB_v1.accdb"

# pass MS Access file path to connection string

db <- odbcConnectAccess2007(file_path)

TABS <- sqlTables(db)$TABLE_NAME

df <- sqlFetch(db, "Query_Survival_Analysis_CompositeCalls")



################# SET FLAG #####################
df$flag <- ifelse(df$PROTOCOLLO != "BO2005" & df$PROTOCOLLO != "EMN02" & 
                    df$TTP_I_months <= 24 & df$Prog_I==1 &
                    df$TP53_adj >= 1.9 & df$MUT_p53_D_SUB_CLON ==0 , 
                  1, 0)
table(df$flag)

table(df$patient_flag, df$flag)

df <- df %>% filter(flag != 1)


#________________ data prep _________________________
# create the label for del tp53 call

df$call50_del_TP53 <- ifelse(df$TP53_adj <= 1.5, 1,0)
df$call40_del_TP53 <- ifelse(df$TP53_adj <= 1.6, 1,0)
df$call30_del_TP53 <- ifelse(df$TP53_adj <= 1.7, 1,0)
df$call20_del_TP53 <- ifelse(df$TP53_adj <= 1.8, 1,0)
df$call10_del_TP53 <- ifelse(df$TP53_adj <= 1.9, 1,0)



#______ revision protocol ________
df$PROTOCOL <- df$PROTOCOLLO %>% recode("FORTE"="AMB", "FP - EMN02"="AMB", "FP - BO2005"="AMB" )
df$PROTOCOL %>% table


#________ compute PFS2 ________ 
# PFS2 is defined as thetime from randomisation (or registration, in non-randomised trials) 
# to second objective disease progression, or death from any cause, whichever first. 
# https://www.ema.europa.eu/en/documents/scientific-guideline/appendix-1-guideline-evaluation-anticancer-medicinal-products-man-methodological-consideration-using_en.pdf

# fill missing ProgII

df$Prog_II[is.na(df$Prog_II)] <- 0


# calcolo nuova variabile = secondPFS_event
df$secondPFS_event <- ifelse( df$Prog_II == 1 | df$OS_event_death == 1, 1, 0)
df$secondPFS_event %>% table

# convert the dates in correct Date format
df$Date0_of_start_IND_therapy <- as.Date(df$Date0_of_start_IND_therapy, origin="1970-01-01")
df$data_II_relapse <- as.Date(df$data_II_relapse, origin="1970-01-01")
df$Date_of_death <- as.Date(df$Date_of_death, origin="1970-01-01")
df$LAST_FULLOW_UP <- as.Date(df$LAST_FULLOW_UP, origin="1970-01-01")
df$PFS_I_date <- as.Date(df$PFS_I_date, origin="1970-01-01")

# calcolo nuova variabile = secondPFS_date
df$secondPFS_date <- dplyr::if_else(df$Prog_II == 1, 
                                    as.Date(df$data_II_relapse, origin="1970-01-01"), 
                                    dplyr::if_else(df$OS_event_death == 1,
                                                   as.Date(df$Date_of_death, origin="1970-01-01"), 
                                                   as.Date(df$LAST_FULLOW_UP,origin="1970-01-01")
                                    ))

df$twoPFS_event <- ifelse(df$Prog_II==1,1, ifelse(df$OS_event_death==1,1,0))


# calulate the twoPFS_time
df$twoPFS_time <- df$secondPFS_date - df$Date0_of_start_IND_therapy

# calulate the twoPFS_months
df$twoPFS_months <- round(df$twoPFS_time / 30.5, 0) %>% as.numeric

df %>% select(LAST_FULLOW_UP,
              Date0_of_start_IND_therapy, 
              PFS_I_date, PFS_I_months, data_II_relapse,
              OS_MESI, secondPFS_date, 
              twoPFS_time, 
              twoPFS_months, 
              Date_of_death, 
              twoPFS_event, 
              OS_event_death, Prog_I,
              Prog_II) %>% View




#_____ cut survival _______
cutoff <- 100
# 
# df$OS_event_death[df$OS_MESI>cutoff] <- 0
# df$OS_MESI[df$OS_MESI>cutoff] <- cutoff
# 
# df$PFS_I_event[df$PFS_I_months>cutoff] <- 0
# df$PFS_I_months[df$PFS_I_months>cutoff] <- cutoff
# 
df$twoPFS_event[df$twoPFS_months>cutoff] <- 0
df$twoPFS_months[df$twoPFS_months>cutoff] <- cutoff



#_____ create *NEW* PFS and OS time variables ________

cut <- 100

df$PFS_date <- dplyr::if_else(df$Prog_I == 1, 
                              as.Date(df$DATA_I_Prog, origin="1970-01-01"), 
                              dplyr::if_else(df$OS_event_death == 1,
                                             as.Date(df$Date_of_death, origin="1970-01-01"), 
                                             as.Date(df$LAST_FULLOW_UP,origin="1970-01-01")
                              ))

df$PFS_time <- df$PFS_date - df$Date0_of_start_IND_therapy

df$PFS_time_months <- df$PFS_time/30.5 %>% as.numeric

df$PFS_I_event[df$PFS_time_months>cut] <- 0 # CUT event 
df$PFS_time_months[df$PFS_time_months>cut] <- cut # CUT 95 months






df$OS_date <- dplyr::if_else(df$OS_event_death == 1, 
                             as.Date(df$Date_of_death, origin="1970-01-01"), as.Date(df$LAST_FULLOW_UP,origin="1970-01-01"))

df$OS_time <- df$OS_date - df$Date0_of_start_IND_therapy

df$OS_time_months <- df$OS_time/30.5 %>% as.numeric

df$OS_event_death[df$OS_time_months>cut] <- 0 # CUT event 
df$OS_time_months[df$OS_time_months>cut] <- cut # CUT 95 months




# df %>% select(PFS_time, PFS_time_months, PFS_I_months, OS_time, OS_time_months, OS_MESI) %>% View




# OS <- Surv( time = df$OS_MESI, event = df$OS_event_death)
# PFS <- Surv( time = df$PFS_I_months, event = df$PFS_I_event)
# twoPFS <- Surv( time=df$twoPFS_months, event = df$twoPFS_event)


OS <- Surv( time = df$OS_time_months, event = df$OS_event_death)
PFS <- Surv( time = df$PFS_time_months, event = df$PFS_I_event)
twoPFS <- Surv( time=df$twoPFS_months, event = df$twoPFS_event)



#____ create analysis variables______
df$TP53_mutation_no_del <- ifelse(df$call10_del_TP53==0 & df$MUT_p53_D_SUB_CLON==1, 1, 0)
df$TP53_del_no_mut <- ifelse(df$call10_del_TP53==1 & df$MUT_p53_D_SUB_CLON==0, 1, 0)
df$del_only_AND_mut_only <- ifelse(df$call10_del_TP53==1 & df$MUT_p53_D_SUB_CLON==0 | df$call10_del_TP53==0 & df$MUT_p53_D_SUB_CLON==1, 1, 0 )
df$Double_Hit <- ifelse(df$MUT_p53_D_SUB_CLON==1 & df$call10_del_TP53==1, 1, 0)

df$group <- ifelse(df$call10_del_TP53==1 & df$MUT_p53_D_SUB_CLON==1, "Double_Hit", 
                   ifelse(df$call10_del_TP53==0 & df$MUT_p53_D_SUB_CLON==0, "WT", "One_Hit")) %>% as.factor



df$group_clonal_subclona <- ifelse(df$TP53_adj<1.5,"del clonal", ifelse(df$TP53_adj<1.9,"subclonal", "wt"))

df$group_clonal_subclona %>% table


df <- df %>% rename(TX = TX_I_line_yes_1_no_0)


#########################################################################
# set outpath
outpath <- "C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/BCJ_REVISIONduefor140821/"

# univariate TX

df_M <- df %>% filter(!is.na(ISS))

OS_M <- Surv( time = df_M$OS_time_months, event = df_M$OS_event_death)
PFS_M <- Surv( time = df_M$PFS_time_months, event = df_M$PFS_I_event)

#____________PFS_______________

write_tsv(data.frame("univariate"),paste0(outpath,"report_cox_TX_290621.txt"), append = T)
with(df_M, coxph(PFS_M ~ TX + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_TX_290621.txt"), append = T, col_names = T)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_TX_290621.txt"), append = T)

write_tsv(data.frame("multivariate"),paste0(outpath,"report_cox_TX_290621.txt"), append = T)
with(df_M, coxph(PFS_M ~ TX + Double_Hit + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .) %>% write_tsv(paste0(outpath,"report_cox_TX_290621.txt"), append = T, col_names = T)


#_______________ OS _______________
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_TX_290621.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_TX_290621.txt"), append = T)

write_tsv(data.frame("univariate"),paste0(outpath,"report_cox_TX_290621.txt"), append = T)
with(df_M, coxph(OS_M ~ TX + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_TX_290621.txt"), append = T, col_names = T)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_TX_290621.txt"), append = T)

write_tsv(data.frame("multivariate"),paste0(outpath,"report_cox_TX_290621.txt"), append = T)
with(df_M, coxph(OS_M ~ TX + Double_Hit + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .) %>% write_tsv(paste0(outpath,"report_cox_TX_290621.txt"), append = T, col_names = T)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_TX_290621.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_TX_290621.txt"), append = T)


# MY MV PFS
with(df_M, coxph(PFS_M ~ TX + Double_Hit + Del_17p_composite + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_multivariate_170221.txt"), append = T, col_names = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_multivariate_170221.txt"), append = T)
with(df_M, coxph(PFS_M ~ TX + Double_Hit + Del_17p_composite + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_multivariate_170221.txt"), append = T, col_names = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_multivariate_170221.txt"), append = T)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_TX_290621.txt"), append = T)

# MY MV OS
write_tsv(data.frame("\n"),paste0(outpath,"report_multivariate_170221.txt"), append = T)
with(df_M, coxph(OS_M ~ TX + plt_m_150 + Double_Hit + Del_17p_composite + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .) %>% write_tsv(paste0(outpath,"report_multivariate_170221.txt"), append = T, col_names = T)

