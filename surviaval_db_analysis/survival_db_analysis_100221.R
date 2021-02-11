
library(RODBC)
library(data.table)
library(tidyverse)
library(survival)
library(survminer)

# full file path to Access DB
file_path <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/TP53_DB_v1.accdb"

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



#_____ cut survival _______

cutoff <- 96

df$OS_event_death[df$OS_MESI>cutoff] <- 0
df$OS_MESI[df$OS_MESI>96] <- cutoff

df$PFS_I_event[df$PFS_I_months>cutoff] <- 0
df$PFS_I_months[df$PFS_I_months>cutoff] <- cutoff


OS <- Surv( time = df$OS_MESI, event = df$OS_event_death)
PFS <- Surv( time = df$PFS_I_months, event = df$PFS_I_event)


##################################################################################
######################### NEW SURV CURVES for PAPER ##############################
##################################################################################

outpath <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/Paper_figures/after_revision_1_110221/"


gmodels::CrossTable(df$call10_del_TP53, df$MUT_p53_D_SUB_CLON, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)


#_____________ 1a) 97 pts (no del - no mut) vs 5 (only mut) _____________

df$TP53_mutation_no_del <- ifelse(df$call10_del_TP53==0 & df$MUT_p53_D_SUB_CLON==1, 1, 0)

df2 <- df %>% filter(call10_del_TP53 == 0) # only no del - 97 (no del - no mut) vs 5 (only mut)
OS2 <- Surv( time = df2$OS_MESI, event = df2$OS_event_death)
PFS2 <- Surv( time = df2$PFS_I_months, event = df2$PFS_I_event)

gg <- ggsurvplot(survfit(OS2 ~ df2$TP53_mutation_no_del, data = df2), pval = T, risk.table = T, xlab = "OS",
                 surv.median.line = "hv", break.time.by = 12, legend= c(0.85,0.9), legend.title="", 
                 legend.labs=c("TP53 wt", "TP53 1 hit (mut+ del-)"), tables.y.text = F, 
                 risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "only_mut_OS.png", path = outpath, dpi = 300, height = 6, width = 8)


gg <- ggsurvplot(survfit(PFS2 ~ df2$TP53_mutation_no_del, data = df2), pval = T, risk.table = T, xlab = "PFS",
                 surv.median.line = "hv", break.time.by = 12, legend= c(0.85,0.9), legend.title="", 
                 legend.labs=c("TP53 wt", "TP53 1 hit (mut+ del-)"), tables.y.text = F, 
                 risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "only_mut_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)





#_____________ 1b) 97 pts (no del - no mut) vs 39 (only del) _____________

df$TP53_del_no_mut <- ifelse(df$call10_del_TP53==1 & df$MUT_p53_D_SUB_CLON==0, 1, 0)

df3 <- df %>% filter(MUT_p53_D_SUB_CLON == 0) 
OS3 <- Surv( time = df3$OS_MESI, event = df3$OS_event_death)
PFS3 <- Surv( time = df3$PFS_I_months, event = df3$PFS_I_event)

gg <- ggsurvplot(survfit(OS3 ~ df3$TP53_del_no_mut, data = df3), pval = T, risk.table = T, xlab = "OS",
                 surv.median.line = "hv", break.time.by = 12, legend= c(0.85,0.9), legend.title="", 
                 legend.labs=c("TP53 wt", "TP53 1 hit (mut- del+)"), tables.y.text = F, 
                 risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "only_del_OS.png", path = outpath, dpi = 300, height = 6, width = 8)


gg <- ggsurvplot(survfit(PFS3 ~ df3$TP53_del_no_mut, data = df3), pval = T, risk.table = T, xlab = "PFS",
                 surv.median.line = "hv", break.time.by = 12, legend= c(0.85,0.9), legend.title="", 
                 legend.labs=c("TP53 wt", "TP53 1 hit (mut- del+)"), tables.y.text = F, 
                 risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "only_del_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)






# 2) Single HIT ---- only mut + only del (39 pts) vs no del no mut (97 pts)

df$del_only_AND_mut_only <- ifelse(df$call10_del_TP53==1 & df$MUT_p53_D_SUB_CLON==0 | df$call10_del_TP53==0 & df$MUT_p53_D_SUB_CLON==1, 1, 0 )

df4 <- df %>% filter(!(MUT_p53_D_SUB_CLON==1 & call10_del_TP53==1)) # only 1 hit - 97 (no del and no mut) vs 39 (only del or only mut)
gmodels::CrossTable(df4$call10_del_TP53, df4$MUT_p53_D_SUB_CLON, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)


OS4 <- Surv( time = df4$OS_MESI, event = df4$OS_event_death)
PFS4 <- Surv( time = df4$PFS_I_months, event = df4$PFS_I_event)

gg <- ggsurvplot(survfit(OS4 ~ df4$del_only_AND_mut_only, data = df4), legend.labs=c("TP53 wt", "TP53 1 hit (mut+ or del+)"),
           pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.85,0.9), 
           legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))

print(gg)

ggsave(plot = print(gg), filename = "1HIT_OS.png", path = outpath, dpi = 300, height = 6, width = 8)


gg <- ggsurvplot(survfit(PFS4 ~ df4$del_only_AND_mut_only, data = df4), legend.labs=c("TP53 wt", "TP53 1 hit (mut+ or del+)"),
           pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.85,0.9), 
           legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))

print(gg)

ggsave(plot = print(gg), filename = "1HIT_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)




# 3) double hit (7 pts) vs no del no mut (97 pts)

df$Double_Hit <- ifelse(df$MUT_p53_D_SUB_CLON==1 & df$call10_del_TP53==1, 1, 0)


df5 <- df %>% filter(!(MUT_p53_D_SUB_CLON==1 & call10_del_TP53==0 | MUT_p53_D_SUB_CLON==0 & call10_del_TP53==1)) # only no del - 97 (no del - no mut) vs 5 (only mut)

OS5 <- Surv( time = df5$OS_MESI, event = df5$OS_event_death)
PFS5 <- Surv( time = df5$PFS_I_months, event = df5$PFS_I_event)

df5$Double_Hit <- ifelse(df5$MUT_p53_D_SUB_CLON==1 & df5$call10_del_TP53==1, 1, 0)

gg <- ggsurvplot(survfit(OS5 ~ df5$Double_Hit, data = df5), legend.labs=c("TP53 wt", "TP53 double-hit"),
           pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.85,0.9), 
           legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))

print(gg)

ggsave(plot = print(gg), filename = "Double-HIT_OS.png", path = outpath, dpi = 300, height = 6, width = 8)


gg <-ggsurvplot(survfit(PFS5 ~ df5$Double_Hit, data = df5), legend.labs=c("TP53 wt", "TP53 double-hit"),
           pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.85,0.9), 
           legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))

print(gg)

ggsave(plot = print(gg), filename = "Double-HIT_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)



# 4) three curves

df$group <- ifelse(df$call10_del_TP53==1 & df$MUT_p53_D_SUB_CLON==1, "Double_Hit", 
                   ifelse(df$call10_del_TP53==0 & df$MUT_p53_D_SUB_CLON==0, "WT", "One_Hit")) %>% as.factor


gg <- ggsurvplot(survfit(OS ~ df$group, data = df) , pval = T, risk.table = T, xlab = "OS",
           surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", 
           legend.labs=c("TP53 double-hit", "TP53 1 Hit", "TP53 wt"), tables.y.text = F, 
           risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "three_curves_OS.png", path = outpath, dpi = 300, height = 6, width = 8)


gg <- ggsurvplot(survfit(PFS ~ df$group, data = df), pval = T, risk.table = T, xlab = "PFS",
           surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", 
           legend.labs=c("TP53 double-hit", "TP53 1 Hit", "TP53 wt"), tables.y.text = F, 
           risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "three_curves_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


#=============== protocollo ===================



gg <- ggsurvplot(survfit(OS ~ df$PROTOCOL, data = df) , pval = T, risk.table = T, xlab = "OS",
                 surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", 
                 tables.y.text = F, 
                 risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "PROTOCOL_OS.png", path = outpath, dpi = 300, height = 6, width = 8)


gg <- ggsurvplot(survfit(PFS ~ df$PROTOCOL, data = df), pval = T, risk.table = T, xlab = "PFS",
                 surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", 
                 tables.y.text = F, 
                 risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "PROTOCOL_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


df$group = relevel(df$group, ref = "WT")
coxph(OS ~ df$group + strata(PROTOCOL), data = df)  %>% summary
coxph(PFS ~ df$group + strata(ISS), data = df) %>% summary


##################### UNIVARIATE ########################


# only mut
coxph(OS2 ~ df2$TP53_mutation_no_del + strata(ISS) + strata(PROTOCOL), data = df2) %>% summary
coxph(PFS2 ~ df2$TP53_mutation_no_del + strata(ISS) + strata(PROTOCOL), data = df2) %>% summary

# only del
coxph(OS3 ~ df3$TP53_del_no_mut + strata(ISS) + strata(PROTOCOL), data = df3) %>% summary
coxph(PFS3 ~ df3$TP53_del_no_mut + strata(ISS) + strata(PROTOCOL), data = df3) %>% summary

# single hit
coxph(OS4 ~ df4$del_only_AND_mut_only + strata(ISS) + strata(PROTOCOL), data = df4) %>% summary
coxph(PFS4 ~ df4$del_only_AND_mut_only + strata(ISS) + strata(PROTOCOL), data = df4) %>% summary

# Double hit
coxph(OS5 ~ df5$Double_Hit + strata(ISS) + strata(PROTOCOL), data = df5) %>% summary
coxph(PFS5 ~ df5$Double_Hit + strata(ISS) + strata(PROTOCOL), data = df5) %>% summary



