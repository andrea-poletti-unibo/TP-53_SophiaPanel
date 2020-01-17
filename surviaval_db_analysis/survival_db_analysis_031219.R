library(RODBC) 
library(data.table)
library(tidyverse)
library(survival)
library(survminer)


db <- odbcConnectAccess2007("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/TP53_DB_v1.accdb")
# sqlTables(db)
df <- sqlFetch(db, "Query_Survival_Analysis")




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

# =======
