# library(RODBC) 
library(odbc)
library(DBI)
library(data.table)
library(tidyverse)
library(survival)
library(survminer)

# devtools::install_version('odbc', '1.2.2', repos="https://cran.rstudio.com/" )

# db <- odbcConnectAccess2007("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/TP53_DB_v1.accdb")
# sqlTables(db)
# df <- sqlFetch(db, "Query_Survival_Analysis")


odbcListDrivers()

# full file path to Access DB
file_path <- "C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/TP53_DB_v1.accdb"

# pass MS Access file path to connection string
accdb_con <- dbConnect(drv = odbc(), .connection_string = paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=",file_path,";"))

dbListTables(accdb_con)

df <- DBI::dbReadTable(accdb_con, "Query_Survival_Analysis")


df$flag <- ifelse(df$PROTOCOLLO != "BO2005" & df$PROTOCOLLO != "EMN02" & 
                    df$TTP_I_months <= 24 & df$Prog_I==1 &
                    df$TP53_adj >= 1.9 & df$MUT_p53_D_SUB_CLON ==0 , 
                  1, 0)

table(df$flag)

df <- df %>% filter(flag != 1)



#________________ data prep _________________________

# create the label for del tp53 call

df$call50_del_TP53 <- ifelse(df$TP53_adj <= 1.5, 1,0)
table(df$call50_del_TP53, df$Fish_DEL_17)

df$call40_del_TP53 <- ifelse(df$TP53_adj <= 1.6, 1,0)
table(df$call40_del_TP53, df$Fish_DEL_17)

df$call30_del_TP53 <- ifelse(df$TP53_adj <= 1.7, 1,0)
table(df$call30_del_TP53, df$Fish_DEL_17)

df$call20_del_TP53 <- ifelse(df$TP53_adj <= 1.8, 1,0)
table(df$call20_del_TP53, df$Fish_DEL_17)

df$call10_del_TP53 <- ifelse(df$TP53_adj <= 1.9, 1,0)
table(df$call10_del_TP53, df$Fish_DEL_17)

#============ Survival analysis =====================

OS <- Surv( time = df$OS_MESI, event = df$OS_event_death)
PFS <- Surv( time = df$PFS_I_months, event = df$PFS_I_event)

#______ t(4;14) ________

ggsurvplot(survfit(OS ~ df$Fish_T_4_14, data = df), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ df$Fish_T_4_14, data = df), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ df$Fish_T_4_14, data = df) %>% summary
coxph(PFS ~ df$Fish_T_4_14, data = df) %>% summary


#______ call50_del_TP53 ________

ggsurvplot(survfit(OS ~ df$call50_del_TP53, data = df), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ df$call50_del_TP53, data = df), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ df$call50_del_TP53 + strata(ISS) , data = df) %>% summary
coxph(PFS ~ df$call50_del_TP53 + strata(ISS), data = df) %>% summary


#______ call40_del_TP53 ________

ggsurvplot(survfit(OS ~ df$call40_del_TP53, data = df), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ df$call40_del_TP53, data = df), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ df$call40_del_TP53+ strata(ISS), data = df) %>% summary
coxph(PFS ~ df$call40_del_TP53+ strata(ISS), data = df) %>% summary

#______ call30_del_TP53 ________

ggsurvplot(survfit(OS ~ df$call30_del_TP53, data = df), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ df$call30_del_TP53, data = df), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ df$call30_del_TP53+ strata(ISS), data = df) %>% summary
coxph(PFS ~ df$call30_del_TP53+ strata(ISS), data = df) %>% summary

#______ call20_del_TP53 ________

ggsurvplot(survfit(OS ~ df$call20_del_TP53, data = df), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ df$call20_del_TP53, data = df), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ df$call20_del_TP53+ strata(ISS), data = df) %>% summary
coxph(PFS ~ df$call20_del_TP53+ strata(ISS), data = df) %>% summary

#______ call10_del_TP53 ________

ggsurvplot(survfit(OS ~ df$call10_del_TP53, data = df), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ df$call10_del_TP53, data = df), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ df$call10_del_TP53+ strata(ISS), data = df) %>% summary
coxph(PFS ~ df$call10_del_TP53+ strata(ISS), data = df) %>% summary


#============== MUT DEL ======================

mutdel <- df %>% filter( TP53_adj < 1.9 & MUT_p53_D_SUB_CLON == 1)

df$mut_del_event <- ifelse(df$TP53_adj < 1.9 & df$MUT_p53_D_SUB_CLON == 1, 1, 0)
table(df$mut_del_event)


homdel <- df %>% filter( TP53_adj < 0.6 )


df$homdel <- ifelse(df$TP53_adj < 0.6, 1, 0)
table(df$homdel)

#______ mutdel survival ________

ggsurvplot(survfit(OS ~ df$mut_del_event, data = df), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ df$mut_del_event, data = df), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ df$mut_del_event, data = df) %>% summary
coxph(PFS ~ df$mut_del_event, data = df) %>% summary


#______ homdel survival ________

ggsurvplot(survfit(OS ~ df$homdel, data = df), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ df$homdel, data = df), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ df$homdel, data = df) %>% summary
coxph(PFS ~ df$homdel, data = df) %>% summary


#______ double inactivation on TP53 _________

df$doubleHIT <- ifelse(df$mut_del_event ==1 | df$homdel == 1, 1,0)
df$doubleHIT %>% table

ggsurvplot(survfit(OS ~ df$doubleHIT, data = df), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ df$doubleHIT, data = df), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ df$doubleHIT, data = df) %>% summary
coxph(PFS ~ df$doubleHIT, data = df) %>% summary


table(df$mut_del_event, df$homdel)


#______ double inactivation with LOH on TP53 _________

df$mut_LOH <- ifelse(df$LOH_call == 1 & df$MUT_p53_D_SUB_CLON==1, 1,0)
df$double_del_LOH <- ifelse(df$mut_del_event ==1 | df$mut_LOH==1, 1,0)

table(df$double_del_LOH, df$mut_LOH)

(df$LOH_call == 1 & df$MUT_p53_D_SUB_CLON==1) %>% table


ggsurvplot(survfit(OS ~ df$double_del_LOH, data = df), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ df$double_del_LOH, data = df), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ df$double_del_LOH, data = df) %>% summary
coxph(PFS ~ df$double_del_LOH, data = df) %>% summary

table(df$doubleHIT, df$double_del_LOH)



#_______________


ggsurvplot(survfit(OS ~ df$MUT_p53_D_SUB_CLON, data = df), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ df$MUT_p53_D_SUB_CLON, data = df), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ df$MUT_p53_D_SUB_CLON, data = df) %>% summary
coxph(PFS ~ df$MUT_p53_D_SUB_CLON, data = df) %>% summary


table(df$doubleHIT, df$double_del_LOH)

GG <- df %>% filter(MUT_p53_D_SUB_CLON==1)

GG$PROTOCOLLO %>% table
mutdel$PROTOCOLLO %>% table

GG %>% select(PROTOCOLLO, PFS_I_months, OS_MESI) %>% View

#================

# CAVO request 09/02/2021 

# 1) 97 pts (no del - no mut) vs 5 (only mut)

df$TP53_mutation_no_del <- ifelse(df$call10_del_TP53==0 & df$MUT_p53_D_SUB_CLON==1, 1,0)

gmodels::CrossTable(df$call10_del_TP53, df$MUT_p53_D_SUB_CLON, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)

df$call10_del_TP53 %>% sum


df3 <- df %>% filter(call10_del_TP53 == 0) # only no del - 97 (no del - no mut) vs 5 (only mut)


OS3 <- Surv( time = df3$OS_MESI, event = df3$OS_event_death)
PFS3 <- Surv( time = df3$PFS_I_months, event = df3$PFS_I_event)


ggsurvplot(survfit(OS3 ~ df3$TP53_mutation_no_del, data = df3), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS3 ~ df3$TP53_mutation_no_del, data = df3), pval = T, risk.table = T, xlab = "PFS")


# 2) only mut + only del (39 pts) vs no del no mut (97 pts)

df4 <- df %>% filter(!(MUT_p53_D_SUB_CLON==1 & call10_del_TP53==1)) # only no del - 97 (no del - no mut) vs 5 (only mut)
gmodels::CrossTable(df4$call10_del_TP53, df4$MUT_p53_D_SUB_CLON, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)


OS4 <- Surv( time = df4$OS_MESI, event = df4$OS_event_death)
PFS4 <- Surv( time = df4$PFS_I_months, event = df4$PFS_I_event)


df4$del_only_AND_mut_only <- ifelse(df4$call10_del_TP53==1 & df4$MUT_p53_D_SUB_CLON==0 |df4$call10_del_TP53==0 & df4$MUT_p53_D_SUB_CLON==1, 1, 0 )

df4$del_only_AND_mut_only %>% sum


ggsurvplot(survfit(OS4 ~ df4$del_only_AND_mut_only, data = df4), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS4 ~ df4$del_only_AND_mut_only, data = df4), pval = T, risk.table = T, xlab = "PFS")



# 3) double hit (7 pts) vs no del no mut (97 pts)

df5 <- df %>% filter(!(MUT_p53_D_SUB_CLON==1 & call10_del_TP53==0 | MUT_p53_D_SUB_CLON==0 & call10_del_TP53==1)) # only no del - 97 (no del - no mut) vs 5 (only mut)
gmodels::CrossTable(df5$call10_del_TP53, df5$MUT_p53_D_SUB_CLON, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)


OS5 <- Surv( time = df5$OS_MESI, event = df5$OS_event_death)
PFS5 <- Surv( time = df5$PFS_I_months, event = df5$PFS_I_event)

df5$Double_Hit <- ifelse(df5$MUT_p53_D_SUB_CLON==1 & df5$call10_del_TP53==1, 1, 0)
df5$Double_Hit %>% sum

ggsurvplot(survfit(OS5 ~ df5$Double_Hit, data = df5), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS5 ~ df5$Double_Hit, data = df5), pval = T, risk.table = T, xlab = "PFS")


# 4) three curves

df$group <- ifelse(df$call10_del_TP53==1 & df$MUT_p53_D_SUB_CLON==1, "Double_Hit", ifelse(df$call10_del_TP53==0 & df$MUT_p53_D_SUB_CLON==0, "WT", "One_Hit"))
df$group %>% table



ggsurvplot(survfit(OS ~ df$group, data = df), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ df$group, data = df), pval = T, risk.table = T, xlab = "PFS")
