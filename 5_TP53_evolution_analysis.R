
library(RODBC)
library(data.table)
library(tidyverse)

db <- odbcConnectAccess2007("C:/Users/andre//Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/TP53_DB_v1.accdb")

TABS <- sqlTables(db)$TABLE_NAME
TABS

# query per tutti i pazienti del DB che hanno i dati molecolari di relapse (65)
query_relapse <- sqlFetch(db, "Copia di Query_Survival_Analysis")

# filtro per i pazienti che sono progrediti (2 eccezioni)
prog_true <- query_relapse %>% filter(Prog_I==1)

# completamento della colonna prog_II: inserimento di 0 nei valori mancanti (pazienti che non sono progrediti ma sono morti o ancora vivi, mancava il dato)
prog_true$Prog_II[is.na(prog_true$Prog_II)] <- 0


prog_true$flag <- ifelse(prog_true$PROTOCOLLO != "BO2005" & prog_true$PROTOCOLLO != "EMN02" & 
                           prog_true$TTP_I_months <= 24 & prog_true$Prog_I==1 &
                           prog_true$TP53_adj_D >= 1.9 & prog_true$MUT_p53_D_SUB_CLON ==0 , 
                         1, 0)

prog_true$flag %>% table

DF2 <- prog_true %>% filter(flag==0)

# del 10%
DF2$call10_del_TP53_D <- ifelse(DF2$TP53_adj_D <= 1.9 , 1, 0)
DF2$call10_del_TP53_R <- ifelse(DF2$TP53_adj_R <= 1.9 , 1, 0)

DF2 %>% names()

evoldf <- DF2 %>% select(UPN, 
                         call10_del_TP53_D, TP53_adj_D,
                         call10_del_TP53_R, TP53_adj_R,
                         MUT_p53_D_SUB_CLON, MUT_P53_D_CCF, 
                         MUT_P53_R_SUB_CLON, MUT_P53_CCF_R)

evoldf <- evoldf[complete.cases(evoldf),]

evoldf %>% summarise(del_D= call10_del_TP53_D %>% sum, 
                     del_R= call10_del_TP53_R %>% sum,
                     mut_D= MUT_p53_D_SUB_CLON %>% sum,
                     mut_R= MUT_P53_R_SUB_CLON %>% sum)

18/53
9/53

gmodels::CrossTable(evoldf$call10_del_TP53_D, evoldf$call10_del_TP53_R, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)

gmodels::CrossTable(evoldf$MUT_p53_D_SUB_CLON, evoldf$MUT_P53_R_SUB_CLON, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)



evoldf$DH_D <- ifelse(evoldf$call10_del_TP53_D ==1 & evoldf$MUT_p53_D_SUB_CLON == 1, 1, 0)
evoldf$DH_D %>% table

evoldf$DH_R <- ifelse(evoldf$call10_del_TP53_R ==1 & evoldf$MUT_P53_R_SUB_CLON == 1, 1, 0)
evoldf$DH_R %>% table

gmodels::CrossTable(evoldf$DH_D, evoldf$DH_R, prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)



##################### DIAGNOSIS #######################

query_diagnosis <- sqlFetch(db, "Query_Survival_Analysis")


query_diagnosis$flag <- ifelse(query_diagnosis$PROTOCOLLO != "BO2005" & query_diagnosis$PROTOCOLLO != "EMN02" & 
                                 query_diagnosis$TTP_I_months <= 24 & query_diagnosis$Prog_I==1 &
                                 query_diagnosis$TP53_adj >= 1.9 & query_diagnosis$MUT_p53_D_SUB_CLON ==0 , 
                         1, 0)

query_diagnosis$flag %>% table

DF1 <- query_diagnosis %>% filter(flag==0)

# del 10%
DF1$call10_del_TP53_D <- ifelse(DF1$TP53_adj <= 1.9 , 1, 0)

DF1$MUT_p53_D_SUB_CLON %>% sum

DF2$MUT_P53_R_SUB_CLON %>% sum




############## proportions D vs R ##############

#_____ dels _____
dels <-   data.frame(Phase= c("Diagnosis","Relapse", "Diagnosis","Relapse"),
                     state= c("del","del","normal","normal"),
                     count = c(sum(DF1$call10_del_TP53_D==1), 
                               sum(evoldf$call10_del_TP53_R==1), 
                               sum(DF1$call10_del_TP53_D==0), 
                               sum(evoldf$call10_del_TP53_R==0))
                     )


dels %>% ggplot(aes(Phase,count, fill=state)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=count), vjust=2, color="white", position = position_stack(), size=3.5) +
  geom_text(x=1, y=154, label="del proportion\n28.6%") +
  geom_text(x=2, y=63, label="del proportion\n34.0%") +
  ylim(0,155)+
  ggtitle("TP53 Deletions proportions ")


41/143
18/53


dels %>% ggplot(aes(Phase, count, fill=state)) + 
  geom_bar(stat="identity", position = position_dodge())+
  geom_text(aes(label=count), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5) +
  ggtitle("TP53 Deletions proportions ")

#_____ muts _____

muts <-   data.frame(Phase= c("Diagnosis","Relapse", "Diagnosis","Relapse"),
                     state= c("mut","mut","normal","normal"),
                     count = c(sum(DF1$MUT_p53_D_SUB_CLON==1), 
                               sum(evoldf$MUT_P53_R_SUB_CLON==1), 
                               sum(DF1$MUT_p53_D_SUB_CLON==0), 
                               sum(evoldf$MUT_P53_R_SUB_CLON==0))
)


muts %>% ggplot(aes(Phase,count, fill=state)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=count), vjust=1.5, color="white", position = position_stack(), size=3.5) +
  geom_text(x=1, y=154, label="mut proportion\n8.3%") +
  geom_text(x=2, y=63, label="mut proportion\n17.0%") +
  ylim(0,155)+
  ggtitle("TP53 Mutation proportions ")

12/143
9/53

muts %>% ggplot(aes(Phase, count, fill=state)) + 
  geom_bar(stat="identity", position = position_dodge())+
  geom_text(aes(label=count), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5) +
  ggtitle("TP53 Deletions proportions ")


#_____ MUTS & DELS _____

DF1$TP53_alteration <- ifelse(DF1$call10_del_TP53_D==1 & DF1$MUT_p53_D_SUB_CLON==1, "DH", 
                              ifelse(DF1$call10_del_TP53_D==1 & DF1$MUT_p53_D_SUB_CLON==0,"del_only", 
                                     ifelse(DF1$call10_del_TP53_D==0 & DF1$MUT_p53_D_SUB_CLON==1, "mut_only","normal"))
                              )
DF1$TP53_alteration %>% table


evoldf$TP53_alteration <- ifelse(evoldf$call10_del_TP53_R==1 & evoldf$MUT_P53_R_SUB_CLON==1, "DH", 
                              ifelse(evoldf$call10_del_TP53_R==1 & evoldf$MUT_P53_R_SUB_CLON==0,"del_only", 
                                     ifelse(evoldf$call10_del_TP53_R==0 & evoldf$MUT_P53_R_SUB_CLON==1, "mut_only","normal"))
)
evoldf$TP53_alteration %>% table


mutdels <-data.frame(Phase= c("Diagnosis","Relapse","Diagnosis","Relapse","Diagnosis","Relapse","Diagnosis","Relapse"),
                     state= c("mut_only","mut_only","normal","normal","del_only","del_only","DH","DH")  %>% factor(levels = c("DH","mut_only","del_only","normal")),
                     count = c(sum(DF1$TP53_alteration=="mut_only"), 
                               sum(evoldf$TP53_alteration=="mut_only"), 
                               sum(DF1$TP53_alteration=="normal"), 
                               sum(evoldf$TP53_alteration=="normal"),
                               sum(DF1$TP53_alteration=="del_only"), 
                               sum(evoldf$TP53_alteration=="del_only"), 
                               sum(DF1$TP53_alteration=="DH"), 
                               sum(evoldf$TP53_alteration=="DH") ),
                     perc = c(sum(DF1$TP53_alteration=="mut_only")/143, 
                               sum(evoldf$TP53_alteration=="mut_only")/53, 
                               sum(DF1$TP53_alteration=="normal")/143, 
                               sum(evoldf$TP53_alteration=="normal")/53,
                               sum(DF1$TP53_alteration=="del_only")/143, 
                               sum(evoldf$TP53_alteration=="del_only")/53, 
                               sum(DF1$TP53_alteration=="DH")/143, 
                               sum(evoldf$TP53_alteration=="DH")/53)
                     )

mutdels$perc <-c( (mutdels$count[c(1,3,5,7)]/143) %>% round(2),
                  (mutdels$count[c(2,4,6,8)]/53) %>% round(2) )

mutdels %>% ggplot(aes(Phase,count, fill=state)) + 
  geom_bar(position = "fill", stat="identity") +
  geom_text(aes(label=count), vjust="top", hjust=1, color="white", position = position_fill(), size=3.5) +
  geom_text(aes(label=(paste0("(",perc %>% round(3) %>% `*`(100),"%)"))), 
            vjust="top",hjust=-0.1, color="white", position = position_fill(), size=3.5) +
  ggtitle("TP53 Alterations proportions and frequency") +
  scale_fill_brewer(palette="RdYlBu") +
  scale_y_continuous(labels = scales::percent) +
  ylab("percentage")





     