
library(odbc)
library(DBI)
library(data.table)
library(tidyverse)
library(survival)
library(survminer)



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
df$call40_del_TP53 <- ifelse(df$TP53_adj <= 1.6, 1,0)
df$call30_del_TP53 <- ifelse(df$TP53_adj <= 1.7, 1,0)
df$call20_del_TP53 <- ifelse(df$TP53_adj <= 1.8, 1,0)
df$call10_del_TP53 <- ifelse(df$TP53_adj <= 1.9, 1,0)

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


# 1) 97 pts (no del - no mut) vs 5 (only mut)

df$TP53_mutation_no_del <- ifelse(df$call10_del_TP53==0 & df$MUT_p53_D_SUB_CLON==1, 1,0)

gmodels::CrossTable(df$call10_del_TP53, df$MUT_p53_D_SUB_CLON, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)


df3 <- df %>% filter(call10_del_TP53 == 0) # only no del - 97 (no del - no mut) vs 5 (only mut)
OS3 <- Surv( time = df3$OS_MESI, event = df3$OS_event_death)
PFS3 <- Surv( time = df3$PFS_I_months, event = df3$PFS_I_event)


ggsurvplot(survfit(OS3 ~ df3$TP53_mutation_no_del, data = df3), pval = T, risk.table = T, xlab = "OS", 
           surv.median.line = "hv", break.time.by = 12, legend.title="", legend.labs=c("TP53 wt", "TP53 1 hit (mut+ del-)")) + theme_survminer(legend=c(0.9,0.9), font.legend=c("bold"))

ggsurvplot(survfit(PFS3 ~ df3$TP53_mutation_no_del, data = df3), pval = T, risk.table = T, xlab = "PFS", 
           surv.median.line = "hv", break.time.by = 12, legend.title="", legend.labs=c("TERTp-wt (all)", "TERT-alt (all)")) + theme_survminer(legend=c(0.9,0.9), font.legend=c("bold"))




# 2) only mut + only del (39 pts) vs no del no mut (97 pts)

df4 <- df %>% filter(!(MUT_p53_D_SUB_CLON==1 & call10_del_TP53==1)) # only no del - 97 (no del - no mut) vs 5 (only mut)
gmodels::CrossTable(df4$call10_del_TP53, df4$MUT_p53_D_SUB_CLON, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)


OS4 <- Surv( time = df4$OS_MESI, event = df4$OS_event_death)
PFS4 <- Surv( time = df4$PFS_I_months, event = df4$PFS_I_event)


df4$del_only_AND_mut_only <- ifelse(df4$call10_del_TP53==1 & df4$MUT_p53_D_SUB_CLON==0 |df4$call10_del_TP53==0 & df4$MUT_p53_D_SUB_CLON==1, 1, 0 )


ggsurvplot(survfit(OS4 ~ df4$del_only_AND_mut_only, data = df4), pval = T, risk.table = T, xlab = "OS",
           surv.median.line = "hv", break.time.by = 12, legend.title="", legend.labs=c("TP53 wt", "TP53 1 hit (mut+ or del+)")) + theme_survminer(legend=c(0.9,0.9), font.legend=c("bold"))

ggsurvplot(survfit(PFS4 ~ df4$del_only_AND_mut_only, data = df4), pval = T, risk.table = T, xlab = "PFS",
           surv.median.line = "hv", break.time.by = 12, legend.title="", legend.labs=c("TP53 wt", "TP53 1 hit (mut+ or del+)")) + theme_survminer(legend=c(0.9,0.9), font.legend=c("bold"))





# 3) double hit (7 pts) vs no del no mut (97 pts)

df5 <- df %>% filter(!(MUT_p53_D_SUB_CLON==1 & call10_del_TP53==0 | MUT_p53_D_SUB_CLON==0 & call10_del_TP53==1)) # only no del - 97 (no del - no mut) vs 5 (only mut)
gmodels::CrossTable(df5$call10_del_TP53, df5$MUT_p53_D_SUB_CLON, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)


OS5 <- Surv( time = df5$OS_MESI, event = df5$OS_event_death)
PFS5 <- Surv( time = df5$PFS_I_months, event = df5$PFS_I_event)

df5$Double_Hit <- ifelse(df5$MUT_p53_D_SUB_CLON==1 & df5$call10_del_TP53==1, 1, 0)
df5$Double_Hit %>% sum

ggsurvplot(survfit(OS5 ~ df5$Double_Hit, data = df5), pval = T, risk.table = T, xlab = "OS",
           surv.median.line = "hv", break.time.by = 12, legend.title="", legend.labs=c("TP53 wt", "TP53 double-hit")) + theme_survminer(legend=c(0.9,0.9), font.legend=c("bold"))

ggsurvplot(survfit(PFS5 ~ df5$Double_Hit, data = df5), pval = T, risk.table = T, xlab = "PFS",
           surv.median.line = "hv", break.time.by = 12, legend.title="", legend.labs=c("TP53 wt", "TP53 double-hit")) + theme_survminer(legend=c(0.9,0.9), font.legend=c("bold"))




# 4) three curves

df$group <- ifelse(df$call10_del_TP53==1 & df$MUT_p53_D_SUB_CLON==1, "Double_Hit", 
                   ifelse(df$call10_del_TP53==0 & df$MUT_p53_D_SUB_CLON==0, "WT", "One_Hit")) %>% as.factor

df$group %>% table



gg <- ggsurvplot(survfit(OS ~ df$group, data = df) , pval = T, risk.table = T, xlab = "OS",
           surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", 
           legend.labs=c("TP53 double-hit", "TP53 1 Hit", "TP53 wt"), tables.y.text = F, 
           risk.table.y.text.col = TRUE, font.legend=c("bold"))
ggsave(plot = print(gg), filename = "three_curves_OS.png", path = "C:/Users/mm_gr/Desktop/", dpi = 300, height = 6, width = 8)


gg <- ggsurvplot(survfit(PFS ~ df$group, data = df), pval = T, risk.table = T, xlab = "PFS",
           surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", 
           legend.labs=c("TP53 double-hit", "TP53 1 Hit", "TP53 wt"), tables.y.text = F, 
           risk.table.y.text.col = TRUE, font.legend=c("bold"))
ggsave(plot = print(gg), filename = "three_curves_PFS.png", path = "C:/Users/mm_gr/Desktop/", dpi = 300, height = 6, width = 8)


df$group = relevel(df$group, ref = "WT")
coxph(OS ~ df$group + strata(ISS), data = df)  %>% summary
coxph(PFS ~ df$group + strata(ISS), data = df) %>% summary




