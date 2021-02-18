library(RODBC)
library(data.table)
library(tidyverse)
library(survival)
library(survminer)

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


#_____ cut survival _______
cutoff <- 96

df$OS_event_death[df$OS_MESI>cutoff] <- 0
df$OS_MESI[df$OS_MESI>96] <- cutoff

df$PFS_I_event[df$PFS_I_months>cutoff] <- 0
df$PFS_I_months[df$PFS_I_months>cutoff] <- cutoff


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


OS <- Surv( time = df$OS_MESI, event = df$OS_event_death)
PFS <- Surv( time = df$PFS_I_months, event = df$PFS_I_event)
twoPFS <- Surv( time=df$twoPFS_months, event = df$twoPFS_event)


#____ create analysis variables______
df$TP53_mutation_no_del <- ifelse(df$call10_del_TP53==0 & df$MUT_p53_D_SUB_CLON==1, 1, 0)
df$TP53_del_no_mut <- ifelse(df$call10_del_TP53==1 & df$MUT_p53_D_SUB_CLON==0, 1, 0)
df$del_only_AND_mut_only <- ifelse(df$call10_del_TP53==1 & df$MUT_p53_D_SUB_CLON==0 | df$call10_del_TP53==0 & df$MUT_p53_D_SUB_CLON==1, 1, 0 )
df$Double_Hit <- ifelse(df$MUT_p53_D_SUB_CLON==1 & df$call10_del_TP53==1, 1, 0)

df$group <- ifelse(df$call10_del_TP53==1 & df$MUT_p53_D_SUB_CLON==1, "Double_Hit", 
                   ifelse(df$call10_del_TP53==0 & df$MUT_p53_D_SUB_CLON==0, "WT", "One_Hit")) %>% as.factor



##################################################################################
######################### NEW SURV CURVES for PAPER ##############################
##################################################################################

outpath <- "C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/Paper_figures/after_revision_1_110221/"


gmodels::CrossTable(df$call10_del_TP53, df$MUT_p53_D_SUB_CLON, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)


#_____________ 1a) 97 pts (no del - no mut) vs 5 (only mut) _____________



df2 <- df %>% filter(call10_del_TP53 == 0) # only no del - 97 (no del - no mut) vs 5 (only mut)

OS2 <- Surv( time = df2$OS_MESI, event = df2$OS_event_death)
PFS2 <- Surv( time = df2$PFS_I_months, event = df2$PFS_I_event)
twoPFS2 <- Surv( time=df2$twoPFS_months, event = df2$twoPFS_event)


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


gg <- ggsurvplot(survfit(twoPFS2 ~ df2$TP53_mutation_no_del, data = df2), pval = T, risk.table = T, xlab = "PFS2",
                 surv.median.line = "hv", break.time.by = 12, legend= c(0.85,0.9), legend.title="", 
                 legend.labs=c("TP53 wt", "TP53 1 hit (mut+ del-)"), tables.y.text = F, 
                 risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "only_mut_PFS2.png", path = outpath, dpi = 300, height = 6, width = 8)



# univariate only mut
c_OS2 <- with(df2, coxph(OS2 ~ TP53_mutation_no_del + strata(ISS) )) %>% summary %>% cbind("OS", .[[7]], .[[8]]) %>% as.data.table(keep.rownames = "var" ) %>% .[,-c(2,4,6,7,9,10)]
fwrite(c_OS2 ,paste0(outpath,"report_univariate_170221.txt"), append = T, sep = "\t")

c_PFS2 <- with(df2, coxph(PFS2 ~ TP53_mutation_no_del + strata(ISS) )) %>% summary %>% cbind("PFS", .[[7]], .[[8]])  %>% as.data.table(keep.rownames = "var" ) %>% .[,-c(2,4,6,7,9,10)]
fwrite(c_PFS2 ,paste0(outpath,"report_univariate_170221.txt"), append = T, sep = "\t")

c_twoPFS2 <- with(df2, coxph(twoPFS2 ~ TP53_mutation_no_del + strata(ISS) )) %>% summary %>% cbind("PFS2", .[[7]], .[[8]]) %>% as.data.table(keep.rownames = "var" ) %>% .[,-c(2,4,6,7,9,10)]
fwrite(c_twoPFS2 ,paste0(outpath,"report_univariate_170221.txt"), append = T, sep = "\t")



#_____________ 1b) 97 pts (no del - no mut) vs 39 (only del) _____________


df3 <- df %>% filter(MUT_p53_D_SUB_CLON == 0) 

OS3 <- Surv( time = df3$OS_MESI, event = df3$OS_event_death)
PFS3 <- Surv( time = df3$PFS_I_months, event = df3$PFS_I_event)
twoPFS3 <- Surv( time=df3$twoPFS_months, event = df3$twoPFS_event)


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


gg <- ggsurvplot(survfit(twoPFS3 ~ df3$TP53_del_no_mut, data = df3), pval = T, risk.table = T, xlab = "PFS2",
                 surv.median.line = "hv", break.time.by = 12, legend= c(0.85,0.9), legend.title="", 
                 legend.labs=c("TP53 wt", "TP53 1 hit (mut- del+)"), tables.y.text = F, 
                 risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "only_del_PFS2.png", path = outpath, dpi = 300, height = 6, width = 8)



# univariate only del
c_OS3 <- with(df3, coxph(OS3 ~ TP53_del_no_mut + strata(ISS) )) %>% summary %>% cbind("OS", .[[7]], .[[8]]) %>% as.data.table(keep.rownames = "var" ) %>% .[,-c(2,4,6,7,9,10)]
fwrite(c_OS3 ,paste0(outpath,"report_univariate_170221.txt"), append = T, sep = "\t")

c_PFS3 <- with(df3, coxph(PFS3 ~ TP53_del_no_mut + strata(ISS) )) %>% summary %>% cbind("PFS", .[[7]], .[[8]])  %>% as.data.table(keep.rownames = "var" ) %>% .[,-c(2,4,6,7,9,10)]
fwrite(c_PFS3 ,paste0(outpath,"report_univariate_170221.txt"), append = T, sep = "\t")

c_twoPFS3 <- with(df3, coxph(twoPFS3 ~ TP53_del_no_mut + strata(ISS) )) %>% summary %>% cbind("PFS2", .[[7]], .[[8]]) %>% as.data.table(keep.rownames = "var" ) %>% .[,-c(2,4,6,7,9,10)]
fwrite(c_twoPFS3 ,paste0(outpath,"report_univariate_170221.txt"), append = T, sep = "\t")




#_______________ 2) Single HIT: only mut + only del (39 pts) vs no del no mut (97 pts) __________________


df4 <- df %>% filter(!(MUT_p53_D_SUB_CLON==1 & call10_del_TP53==1)) # only 1 hit - 97 (no del and no mut) vs 39 (only del or only mut)
gmodels::CrossTable(df4$call10_del_TP53, df4$MUT_p53_D_SUB_CLON, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)


OS4 <- Surv( time = df4$OS_MESI, event = df4$OS_event_death)
PFS4 <- Surv( time = df4$PFS_I_months, event = df4$PFS_I_event)
twoPFS4 <- Surv( time = df4$twoPFS_months, event = df4$twoPFS_event)

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


gg <- ggsurvplot(survfit(twoPFS4 ~ df4$del_only_AND_mut_only, data = df4), legend.labs=c("TP53 wt", "TP53 1 hit (mut+ or del+)"),
                 pval = T, risk.table = T, xlab = "PFS2", surv.median.line = "hv", break.time.by = 12, legend= c(0.85,0.9), 
                 legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))

print(gg)

ggsave(plot = print(gg), filename = "1HIT_PFS2.png", path = outpath, dpi = 300, height = 6, width = 8)


# univariate 1 Hit
c_OS4 <- with(df4, coxph(OS4 ~ del_only_AND_mut_only + strata(ISS) )) %>% summary %>% cbind("OS", .[[7]], .[[8]]) %>% as.data.table(keep.rownames = "var" ) %>% .[,-c(2,4,6,7,9,10)]
fwrite(c_OS4 ,paste0(outpath,"report_univariate_170221.txt"), append = T, sep = "\t")

c_PFS4 <- with(df4, coxph(PFS4 ~ del_only_AND_mut_only + strata(ISS) )) %>% summary %>% cbind("PFS", .[[7]], .[[8]])  %>% as.data.table(keep.rownames = "var" ) %>% .[,-c(2,4,6,7,9,10)]
fwrite(c_PFS4 ,paste0(outpath,"report_univariate_170221.txt"), append = T, sep = "\t")

c_twoPFS4 <- with(df4, coxph(twoPFS4 ~ del_only_AND_mut_only + strata(ISS) )) %>% summary %>% cbind("PFS2", .[[7]], .[[8]]) %>% as.data.table(keep.rownames = "var" ) %>% .[,-c(2,4,6,7,9,10)]
fwrite(c_twoPFS4 ,paste0(outpath,"report_univariate_170221.txt"), append = T, sep = "\t")



#_________________ 3) double hit (7 pts) vs no del no mut (97 pts) _________________________



df5 <- df %>% filter(!(MUT_p53_D_SUB_CLON==1 & call10_del_TP53==0 | MUT_p53_D_SUB_CLON==0 & call10_del_TP53==1)) # only no del - 97 (no del - no mut) vs 5 (only mut)

OS5 <- Surv( time = df5$OS_MESI, event = df5$OS_event_death)
PFS5 <- Surv( time = df5$PFS_I_months, event = df5$PFS_I_event)
twoPFS5 <- Surv( time = df5$twoPFS_months, event = df5$twoPFS_event)


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



gg <-ggsurvplot(survfit(twoPFS5 ~ df5$Double_Hit, data = df5), legend.labs=c("TP53 wt", "TP53 double-hit"),
                pval = T, risk.table = T, xlab = "PFS2", surv.median.line = "hv", break.time.by = 12, legend= c(0.85,0.9), 
                legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))

print(gg)

ggsave(plot = print(gg), filename = "Double-HIT_PFS2.png", path = outpath, dpi = 300, height = 6, width = 8)



# Univariate Double-Hit
c_OS5 <- with(df5, coxph(OS5 ~ Double_Hit + strata(ISS) )) %>% summary %>% cbind("OS", .[[7]], .[[8]]) %>% as.data.table(keep.rownames = "var" ) %>% .[,-c(2,4,6,7,9,10)]
fwrite(c_OS5 ,paste0(outpath,"report_univariate_170221.txt"), append = T, sep = "\t")

c_PFS5 <- with(df5, coxph(PFS5 ~ Double_Hit + strata(ISS) )) %>% summary %>% cbind("PFS", .[[7]], .[[8]])  %>% as.data.table(keep.rownames = "var" ) %>% .[,-c(2,4,6,7,9,10)]
fwrite(c_PFS5 ,paste0(outpath,"report_univariate_170221.txt"), append = T, sep = "\t")

c_twoPFS5 <- with(df5, coxph(twoPFS5 ~ Double_Hit + strata(ISS) )) %>% summary %>% cbind("PFS2", .[[7]], .[[8]]) %>% as.data.table(keep.rownames = "var" ) %>% .[,-c(2,4,6,7,9,10)]
fwrite(c_twoPFS5 ,paste0(outpath,"report_univariate_170221.txt"), append = T, sep = "\t")



#____________________ 4) three curves _______________________

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



gg <- ggsurvplot(survfit(twoPFS ~ df$group, data = df), pval = T, risk.table = T, xlab = "PFS2",
                 surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", 
                 legend.labs=c("TP53 double-hit", "TP53 1 Hit", "TP53 wt"), tables.y.text = F, 
                 risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "three_curves_PFS2.png", path = outpath, dpi = 300, height = 6, width = 8)


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


# df$group = relevel(df$group, ref = "WT")
# coxph(OS ~ df$group + strata(PROTOCOL), data = df)  %>% summary
# coxph(PFS ~ df$group + strata(ISS), data = df) %>% summary


####################### UNIVARIATE ########################
# 
# # only mut
# c_OS2 <- with(df2, coxph(OS2 ~ TP53_mutation_no_del + strata(ISS) + strata(PROTOCOL))) %>% summary %>% cbind(.[[7]], .[[8]]) %>% as.data.frame()
# c_OS2
# 
# c_PFS2 <- with(df2, coxph(PFS2 ~ TP53_mutation_no_del + strata(ISS) + strata(PROTOCOL))) %>% summary %>% cbind(.[[7]], .[[8]]) %>% as.data.frame()
# c_PFS2
# 
# # only del
# coxph(OS3 ~ df3$TP53_del_no_mut + strata(ISS) + strata(PROTOCOL), data = df3) %>% summary
# coxph(PFS3 ~ df3$TP53_del_no_mut + strata(ISS) + strata(PROTOCOL), data = df3) %>% summary
# 
# # single hit
# coxph(OS4 ~ df4$del_only_AND_mut_only + strata(ISS) + strata(PROTOCOL), data = df4) %>% summary
# coxph(PFS4 ~ df4$del_only_AND_mut_only + strata(ISS) + strata(PROTOCOL), data = df4) %>% summary
# 
# # Double hit
# coxph(OS5 ~ df5$Double_Hit + strata(ISS) + strata(PROTOCOL), data = df5) %>% summary
# coxph(PFS5 ~ df5$Double_Hit + strata(ISS) + strata(PROTOCOL), data = df5) %>% summary

coxph(OS ~ df$Del_13q_composite + strata(ISS) , data = df) %>% summary
coxph(PFS ~ df$Del_13q_composite + strata(ISS) , data = df) %>% summary



######################### second PFS analysis ############################


#_____________ compute second PFS _____________


# query per tutti i pazienti del DB che hanno i dati molecolari di relapse (65)
query_relapse <- sqlFetch(db, "Copia di Query_Survival_Analysis")

# filtro per i pazienti che sono progrediti (2 eccezioni)
prog_true <- query_relapse %>% filter(Prog_I==1)

# completamento della colonna prog_II: inserimento di 0 nei valori mancanti (pazienti che non sono progrediti ma sono morti o ancora vivi, mancava il dato)
prog_true$Prog_II[is.na(prog_true$Prog_II)] <- 0

# calcolo nuova variabile = secondPFS_event
prog_true$secondPFS_event <- ifelse( prog_true$Prog_II == 1 | prog_true$OS_event_death == 1, 1, 0)
prog_true$secondPFS_event %>% table

# convert the dates in correct Date format
prog_true$data_II_relapse <- as.Date(prog_true$data_II_relapse, origin="1970-01-01")
prog_true$Date_of_death <- as.Date(prog_true$Date_of_death, origin="1970-01-01")
prog_true$LAST_FULLOW_UP <- as.Date(prog_true$LAST_FULLOW_UP, origin="1970-01-01")
prog_true$PFS_I_date <- as.Date(prog_true$PFS_I_date, origin="1970-01-01")

# calcolo nuova variabile = secondPFS_date
prog_true$secondPFS_date <- dplyr::if_else(prog_true$Prog_II == 1, 
                                       as.Date(prog_true$data_II_relapse, origin="1970-01-01"), 
                                       dplyr::if_else(prog_true$OS_event_death == 1, 
                                                      as.Date(prog_true$Date_of_death, origin="1970-01-01"), 
                                                      as.Date(prog_true$LAST_FULLOW_UP,origin="1970-01-01")
                                       ))

prog_true %>% select(Prog_II, data_II_relapse, OS_event_death, Date_of_death, LAST_FULLOW_UP, secondPFS_event, secondPFS_date) %>% View

# calulate the secondPFS_time
prog_true$secondPFS_time <- prog_true$secondPFS_date - prog_true$PFS_I_date
prog_true %>% select( PFS_I_date, secondPFS_date, secondPFS_time) %>% View

# calulate the secondPFS_months
prog_true$secondPFS_months <- round(prog_true$secondPFS_time / 30.5, 0) 

#______________ set flag and filter________
prog_true$flag <- ifelse(prog_true$PROTOCOLLO != "BO2005" & prog_true$PROTOCOLLO != "EMN02" & 
                           prog_true$TTP_I_months <= 24 & prog_true$Prog_I==1 &
                           prog_true$TP53_adj_D >= 1.9 & prog_true$MUT_p53_D_SUB_CLON ==0 , 
                         1, 0)

DF2 <- prog_true %>% filter(flag==0)


# ISS distribution Check
DF2$ISS %>% table 

# del 10%
DF2$call10_del_TP53_D <- ifelse(DF2$TP53_adj_D <= 1.9 , 1, 0)
DF2$call10_del_TP53_R <- ifelse(DF2$TP53_adj_R <= 1.9 , 1, 0)

DF2$call10_del_TP53_D %>% table
DF2$call10_del_TP53_R %>% table

gmodels::CrossTable(DF2$call10_del_TP53_R, DF2$MUT_P53_R_SUB_CLON, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)

#________________ a) only del (14 pts) vs wt (30 pts) - second PFS __________________

DF2a <- DF2 %>% filter(DF2$MUT_P53_CCF_R==0)

secondPFSa <- Surv( time = DF2a$secondPFS_months, event = DF2a$secondPFS_event)

gg <- ggsurvplot(survfit(secondPFSa ~ DF2a$call10_del_TP53_R, data = DF2a) , pval = T, risk.table = T, xlab = "second PFS",
                 surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", 
                 legend.labs=c("TP53 wt", "TP53 1 hit (mut- del+)"), tables.y.text = F, 
                 risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "only_del_secondPFS.png", path = outpath, dpi = 300, height = 6, width = 8)


# univariata
c_secPFSa <- with(DF2a, coxph(secondPFSa ~ call10_del_TP53_R + strata(ISS) )) %>% summary %>% cbind("secondPFS", .[[7]], .[[8]]) %>% as.data.table(keep.rownames = "var" ) %>% .[,-c(2,4,6,7,9,10)]
fwrite(c_secPFSa ,paste0(outpath,"report_univariate_170221.txt"), append = T, sep = "\t")



#_______________ b) One-Hit (14 pts) vs wt (30 pts) - second PFS ___________________

DF2$del_only_AND_mut_only_R <- ifelse(DF2$call10_del_TP53_R==1 & DF2$MUT_P53_R_SUB_CLON==0 | DF2$call10_del_TP53_R==0 & DF2$MUT_P53_R_SUB_CLON==1, 1, 0 )

# DF2$del_only_AND_mut_only_D <- ifelse(DF2$call10_del_TP53_D==1 & DF2$MUT_p53_D_SUB_CLON==0 | DF2$call10_del_TP53_D==0 & DF2$MUT_p53_D_SUB_CLON==1, 1, 0 )

DF2b <- DF2 %>% filter(DF2$MUT_P53_CCF_R==0 | DF2$call10_del_TP53_R ==0)

gmodels::CrossTable(DF2b$call10_del_TP53_R, DF2b$MUT_P53_R_SUB_CLON, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)
# gmodels::CrossTable(DF2b$call10_del_TP53_D, DF2b$MUT_p53_D_SUB_CLON, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)


secondPFSb <- Surv( time = DF2b$secondPFS_months, event = DF2b$secondPFS_event)


gg <- ggsurvplot(survfit(secondPFSb ~ DF2b$del_only_AND_mut_only_R, data = DF2b) , pval = T, risk.table = T, xlab = "second PFS",
                 surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", 
                 legend.labs=c("TP53 wt", "TP53 1 hit (mut+ or del+)"), tables.y.text = F, 
                 risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "1HIT_secondPFS.png", path = outpath, dpi = 300, height = 6, width = 8)


# univariata
c_secPFSb <- with(DF2b, coxph(secondPFSb ~ del_only_AND_mut_only_R + strata(ISS) )) %>% summary %>% cbind("secondPFS", .[[7]], .[[8]]) %>% as.data.table(keep.rownames = "var" ) %>% .[,-c(2,4,6,7,9,10)]
fwrite(c_secPFSb ,paste0(outpath,"report_univariate_170221.txt"), append = T, sep = "\t")


#_______________ c) Double-Hit (5 pts) vs wt (30 pts) - second PFS ___________________

DF2$Double_Hit_R <- ifelse(DF2$call10_del_TP53_R==1 & DF2$MUT_P53_R_SUB_CLON==1, 1, 0 )

DF2c <- DF2 %>% filter(DF2$MUT_P53_CCF_R==0 & DF2$call10_del_TP53_R ==0 | DF2$Double_Hit_R ==1)

gmodels::CrossTable(DF2c$call10_del_TP53_R, DF2c$MUT_P53_R_SUB_CLON, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)

secondPFSc <- Surv( time = DF2c$secondPFS_months, event = DF2c$secondPFS_event)


gg <- ggsurvplot(survfit(secondPFSc ~ DF2c$Double_Hit_R, data = DF2c) , pval = T, risk.table = T, xlab = "second PFS",
                 surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", 
                 legend.labs=c("TP53 wt", "TP53 double-hit"), tables.y.text = F, 
                 risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "Double-HIT_secondPFS.png", path = outpath, dpi = 300, height = 6, width = 8)


# univariata
c_secPFSc <- with(DF2c, coxph(secondPFSc ~ Double_Hit_R + strata(ISS) )) %>% summary %>% cbind("secondPFS", .[[7]], .[[8]]) %>% as.data.table(keep.rownames = "var" ) %>% .[,-c(2,4,6,7,9,10)]
fwrite(c_secPFSc ,paste0(outpath,"report_univariate_170221.txt"), append = T, sep = "\t")


#_______________ d) Three curves: Double-Hit (5 pts) vs Single-Hit (14 pts) vs wt (30 pts) - second PFS ___________________

DF2$group <- ifelse(DF2$call10_del_TP53_R==1 & DF2$MUT_P53_R_SUB_CLON==1, "Double_Hit", 
                    ifelse(DF2$call10_del_TP53_R==0 & DF2$MUT_P53_R_SUB_CLON==0, "WT", "One_Hit")) %>% as.factor


gmodels::CrossTable(DF2$group, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)

secondPFS <- Surv( time = DF2$secondPFS_months, event = DF2$secondPFS_event)


gg <- ggsurvplot(survfit(secondPFS ~ DF2$group, data = DF2) , pval = T, risk.table = T, xlab = "second PFS",
                 surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", 
                 legend.labs=c("TP53 double-hit", "TP53 1 Hit", "TP53 wt"), tables.y.text = F, 
                 risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "three_curves_secondPFS.png", path = outpath, dpi = 300, height = 6, width = 8)


################################### multivariate ############################################

# EXCLUDE MISSING ISS PTS

df_M <- df %>% filter(!is.na(ISS))

OS_M <- Surv( time = df_M$OS_MESI, event = df_M$OS_event_death)
PFS_M <- Surv( time = df_M$PFS_I_months, event = df_M$PFS_I_event)
# 
# Covariate testate:
# "D_traslocazione","TRASLOCATI","Fish_T_4_14", "Fish_T_6_14_", "Fish_T_11_14",
# "Fish_T_14_16", "Fish_T_14_20",  "Del_13q_composite", "Del_17p_composite", "Del_1p_composite", 
# "Amp_1q_composite", "Fish_iper", "plt_m_150", "Age", "Age_Cat"


#____________PFS_______________

with(df_M, coxph(PFS_M ~ call10_del_TP53 + strata(ISS) )) %>% summary
with(df_M, coxph(PFS_M ~ MUT_p53_D_SUB_CLON + strata(ISS) )) %>% summary
with(df_M, coxph(PFS_M ~ del_only_AND_mut_only + strata(ISS) )) %>% summary
with(df_M, coxph(PFS_M ~ Double_Hit + strata(ISS) )) %>% summary


with(df_M, coxph(PFS_M ~ D_traslocazione + strata(ISS) )) %>% summary
with(df_M, coxph(PFS_M ~ TRASLOCATI + strata(ISS) )) %>% summary
with(df_M, coxph(PFS_M ~ Fish_T_4_14 + strata(ISS) )) %>% summary
with(df_M, coxph(PFS_M ~ Fish_T_6_14_ + strata(ISS) )) %>% summary
with(df_M, coxph(PFS_M ~ Fish_T_11_14 + strata(ISS) )) %>% summary # **
with(df_M, coxph(PFS_M ~ Fish_T_14_16 + strata(ISS) )) %>% summary
with(df_M, coxph(PFS_M ~ Fish_T_14_20 + strata(ISS) )) %>% summary # **
with(df_M, coxph(PFS_M ~ Del_13q_composite + strata(ISS) )) %>% summary # .
with(df_M, coxph(PFS_M ~ Del_17p_composite + strata(ISS) )) %>% summary # **
with(df_M, coxph(PFS_M ~ Del_1p_composite + strata(ISS) )) %>% summary
with(df_M, coxph(PFS_M ~ Amp_1q_composite + strata(ISS) )) %>% summary
with(df_M, coxph(PFS_M ~ Fish_iper + strata(ISS) )) %>% summary
with(df_M, coxph(PFS_M ~ plt_m_150 + strata(ISS) )) %>% summary
with(df_M, coxph(PFS_M ~ Age + strata(ISS) )) %>% summary



with(df_M, coxph(PFS_M ~ call10_del_TP53 + MUT_p53_D_SUB_CLON  + strata(ISS) )) %>% summary
with(df_M, coxph(PFS_M ~ Del_17p_composite + MUT_p53_D_SUB_CLON  + strata(ISS) )) %>% summary
with(df_M, coxph(PFS_M ~ Del_17p_composite + Double_Hit  + strata(ISS) )) %>% summary

# MY MV PFS
with(df_M, coxph(PFS_M ~ Fish_T_11_14 + Fish_T_14_20 + Double_Hit + Del_13q_composite + strata(ISS) )) %>% summary
with(df_M, coxph(PFS_M ~ Fish_T_11_14 + Fish_T_14_20 + Double_Hit + strata(ISS) )) %>% summary


#_______________ OS _______________
del_only_AND_mut_only
with(df_M, coxph(OS_M ~ call10_del_TP53 + strata(ISS) )) %>% summary
with(df_M, coxph(OS_M ~ MUT_p53_D_SUB_CLON + strata(ISS) )) %>% summary
with(df_M, coxph(OS_M ~ del_only_AND_mut_only + strata(ISS) )) %>% summary
with(df_M, coxph(OS_M ~ Double_Hit + strata(ISS) )) %>% summary

with(df_M, coxph(OS_M ~ D_traslocazione + strata(ISS) )) %>% summary
with(df_M, coxph(OS_M ~ TRASLOCATI + strata(ISS) )) %>% summary
with(df_M, coxph(OS_M ~ Fish_T_4_14 + strata(ISS) )) %>% summary
with(df_M, coxph(OS_M ~ Fish_T_6_14_ + strata(ISS) )) %>% summary
with(df_M, coxph(OS_M ~ Fish_T_11_14 + strata(ISS) )) %>% summary # *
with(df_M, coxph(OS_M ~ Fish_T_14_16 + strata(ISS) )) %>% summary
with(df_M, coxph(OS_M ~ Fish_T_14_20 + strata(ISS) )) %>% summary # **
with(df_M, coxph(OS_M ~ Del_13q_composite + strata(ISS) )) %>% summary # ***
with(df_M, coxph(OS_M ~ Del_17p_composite + strata(ISS) )) %>% summary # .
with(df_M, coxph(OS_M ~ Del_1p_composite + strata(ISS) )) %>% summary # .
with(df_M, coxph(OS_M ~ Amp_1q_composite + strata(ISS) )) %>% summary
with(df_M, coxph(OS_M ~ Fish_iper + strata(ISS) )) %>% summary # .
with(df_M, coxph(OS_M ~ plt_m_150 + strata(ISS) )) %>% summary # *
with(df_M, coxph(OS_M ~ Age + strata(ISS) )) %>% summary


with(df_M, coxph(OS_M ~ call10_del_TP53 + Del_13q_composite  + strata(ISS) )) %>% summary
with(df_M, coxph(OS_M ~ Del_17p_composite + plt_m_150 + MUT_p53_D_SUB_CLON  + strata(ISS) )) %>% summary
with(df_M, coxph(OS_M ~ Del_17p_composite + plt_m_150 + Del_13q_composite + Double_Hit  + strata(ISS) )) %>% summary


# MY MV OS

with(df_M, coxph(OS_M ~ Fish_T_11_14 + Fish_T_14_20 + Del_13q_composite + Del_1p_composite + Fish_iper + plt_m_150+ Double_Hit  + strata(ISS) )) %>% summary
with(df_M, coxph(OS_M ~ Fish_T_14_20 + Del_13q_composite + Double_Hit  + strata(ISS) )) %>% summary



########################## REPORT MV ###########################
library(broom)

# del 10% TP53 and mut
MV_PFS_1 <- with(df_M, coxph(PFS_M ~ MUT_p53_D_SUB_CLON + call10_del_TP53  + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)
write_tsv(MV_PFS_1, paste0(outpath,"report_multivariate_170221.txt"), append = T, col_names = T)

MV_OS_1 <- with(df_M, coxph(OS_M ~ MUT_p53_D_SUB_CLON + call10_del_TP53  + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)
write_tsv(MV_PFS_1, paste0(outpath,"report_multivariate_170221.txt"), append = T, col_names = T)


# del 17 and MUT
MV_PFS_2 <- with(df_M, coxph(PFS_M ~ MUT_p53_D_SUB_CLON + Del_17p_composite + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)
write_tsv(MV_PFS_2, paste0(outpath,"report_multivariate_170221.txt"), append = T, col_names = T)

MV_OS_2 <- with(df_M, coxph(OS_M ~ MUT_p53_D_SUB_CLON + Del_17p_composite  + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)
write_tsv(MV_OS_2, paste0(outpath,"report_multivariate_170221.txt"), append = T, col_names = T)



# del 17 and DoubleHit
MV_PFS_3 <- with(df_M, coxph(PFS_M ~ Double_Hit  + Del_17p_composite + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)
write_tsv(MV_PFS_3, paste0(outpath,"report_multivariate_170221.txt"), append = T, col_names = T)

MV_OS_3 <- with(df_M, coxph(OS_M ~ Double_Hit  + Del_17p_composite  + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)
write_tsv(MV_OS_3, paste0(outpath,"report_multivariate_170221.txt"), append = T, col_names = T)



# del 17 and DoubleHit
MV_PFS_3 <- with(df_M, coxph(PFS_M ~ Double_Hit  + Del_17p_composite + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)
write_tsv(MV_PFS_3, paste0(outpath,"report_multivariate_170221.txt"), append = T, col_names = T)

MV_OS_3 <- with(df_M, coxph(OS_M ~ Double_Hit  + Del_17p_composite  + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)
write_tsv(MV_OS_3, paste0(outpath,"report_multivariate_170221.txt"), append = T, col_names = T)




# del 17 - DoubleHit - del 13 - plt m 150
MV_PFS_3 <- with(df, coxph(PFS_M ~ Double_Hit  + Del_17p_composite + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)
write_tsv(MV_PFS_3, paste0(outpath,"report_multivariate_170221.txt"), append = T, col_names = T)

MV_OS_3 <- with(df, coxph(OS_M ~ Double_Hit  +  Del_13q_composite + strata(ISS) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)
write_tsv(MV_OS_3, paste0(outpath,"report_multivariate_170221.txt"), append = T, col_names = T)

