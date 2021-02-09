library(ggplot2)
library(reshape2)
library(dplyr)
library(grid)

#======================= EXAMPLE ==========================
set.seed(33)
dat = data.frame(Subject = 1:10, 
                 Months = sample(4:20, 10, replace=TRUE),
                 Treated = sample(0:1, 10, replace=TRUE),
                 Stage = sample(1:4, 10, replace=TRUE),
                 Continued = sample(0:1, 10, replace=TRUE))

dat = dat %>%
  group_by(Subject) %>%
  mutate(Complete=sample(c(4:(max(Months)-1),NA), 1, 
                         prob=c(rep(1, length(4:(max(Months)-1))),5), replace=TRUE),
         Partial=sample(c(4:(max(Months)-1),NA), 1, 
                        prob=c(rep(1, length(4:(max(Months)-1))),5), replace=TRUE),
         Durable=sample(c(-0.5,NA), 1, replace=TRUE))

# Order Subjects by Months
dat$Subject = factor(dat$Subject, levels=dat$Subject[order(dat$Months)])

# Melt part of data frame for adding points to bars
dat.m = melt(dat %>% select(Subject, Months, Complete, Partial, Durable),
             id.var=c("Subject","Months"))

ggplot(dat, aes(Subject, Months)) +
  geom_bar(stat="identity", aes(fill=factor(Stage)), colour="black", width=0.7) +
  geom_point(data=dat.m, 
             aes(Subject, value, colour=variable, shape=variable), size=4) +
  geom_segment(data=dat %>% filter(Continued==1), 
               aes(x=Subject, xend=Subject, y=Months + 0.1, yend=Months + 1), 
               pch=15, size=0.8, arrow=arrow(type="closed", length=unit(0.1,"in"))) +
  coord_flip() +
  scale_fill_manual(values=hcl(seq(15,375,length.out=5)[1:4],100,70)) +
  scale_colour_manual(values=c(hcl(seq(15,375,length.out=3)[1:2],100,40),"black")) +
  scale_y_continuous(limits=c(-1,20), breaks=0:20) +
  labs(fill="Disease Stage", colour="", shape="", 
       x="Subject Recevied Study Drug") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())


#======================= TP53 SwimmerPlot ==========================

import.pts <- readxl::read_xlsx("../swimmerplot_060919.xlsx", sheet = 2)
import.mut <- readxl::read_xlsx("../swimmerplot_060919.xlsx", sheet = 1)

#____ data prep ____
dat$Subject = factor(dat$Subject, levels=dat$Subject[order(dat$Months)])

data.pts <- import.pts 
data.pts$UPN <- factor(data.pts$UPN, levels = data.pts$UPN[order(data.pts$MESI)]) # REORDER SAMPLES

# filter and join info
data.mut <- import.mut %>%filter(TO_PLOT==1)
data.mut <- left_join(data.mut, data.pts[,-2], by="UPN")

# create a MUT ID
data.mut$MUT_ID <- paste(data.mut$POSITION, data.mut$REF, data.mut$ALT, sep = "_")

# diagnosis MUTs
data.mut_D <- data.mut %>% filter(FASE %in% c("diagnosis", "diagnosis-relapse"))
data.mut_D$CCF <- coalesce(as.numeric(data.mut_D$CCF_D), as.numeric(data.mut_D$VAF_D))
data.mut_D$TIME <- 0
data.mut_D <- data.mut_D %>% select(UPN, FASE,TYPE,GROUP, MUT_ID, PT_MUT_CODE, CCF, TIME)

# relapse MUTs
data.mut_R <- data.mut %>% filter(FASE %in% c("relapse", "diagnosis-relapse"))
data.mut_R$CCF <- coalesce(as.numeric(data.mut_R$CCF_R), as.numeric(data.mut_R$VAF_R))
data.mut_R$TIME <- data.mut_R$MESI
data.mut_R <- data.mut_R %>% select(UPN, FASE,TYPE,GROUP, MUT_ID,CCF, PT_MUT_CODE, TIME)

# merge D and R
plot.mut <- rbind(data.mut_D,data.mut_R)


data.pts %>% ggplot(aes(UPN, MESI)) +
  geom_bar(stat="identity", aes(fill=factor(GROUP)), colour="black", alpha=0.2, width=0.7) +
  geom_segment(data=data.mut %>% filter(Death==0), 
               aes(x=UPN, xend=UPN, y=MESI + 0.5, yend=MESI + 4), 
               pch=15, size=0.8, arrow=arrow(type="closed", length=unit(0.1,"in"))) +
  geom_point(data=plot.mut, 
             aes(UPN, TIME, colour=PT_MUT_CODE, shape=FASE, size=CCF), position=position_dodge(width=1)) +
  coord_flip() +
  scale_fill_manual(values=c("dodgerblue","salmon")) +
  # scale_colour_manual(values=c(hcl(seq(15,375,length.out=3)[1:2],100,40),"black")) +
  labs(title = "TP53 SwimmerPlot", fill="Disease Stage", colour="", shape="", 
       x="Patients with TP53 mutation(s)") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())



############################## new DB swimmerplot 11/03/2020 ##################################

library(RODBC)

db <- odbcConnectAccess2007("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/TP53_DB_v1.accdb")

TABS <- sqlTables(db)$TABLE_NAME
TABS

# query Diagnosis
Q_Dia <- sqlFetch(db, "Query_Survival_Analysis")
# query per tutti i pazienti del DB che hanno i dati molecolari di relapse (65)
Q_Rel <- sqlFetch(db, "Copia di Query_Survival_Analysis")

pts_id <- sqlFetch(db,"valid_pts_ID_table")

#=========== data set up diagnosis ========

Q_Dia$flag <- ifelse(Q_Dia$PROTOCOLLO != "BO2005" & Q_Dia$PROTOCOLLO != "EMN02" & 
                     Q_Dia$TTP_I_months <= 24 & Q_Dia$Prog_I==1 &
                     Q_Dia$TP53_adj >= 1.9 & Q_Dia$MUT_p53_D_SUB_CLON ==0 , 
                     1, 0)

Q_Dia$flag %>% table  

DF1 <- Q_Dia %>% filter(flag==0)

DF1$del_TP53_SNP_D <- ifelse(DF1$TP53_adj<1.9,1,0)

DF1$del_TP53_SNP_D %>% table

#=========== data set up relapse ==========

# filtro per i pazienti che sono progrediti (2 eccezioni)
prog_true <- Q_Rel %>% filter(Prog_I==1)
# completamento della colonna prog_II: inserimento di 0 nei valori mancanti (pazienti che non sono progrediti ma sono morti o ancora vivi, mancava il dato)
prog_true$Prog_II[is.na(prog_true$Prog_II)] <- 0
# calcolo nuova variabile = PFS_2_event
prog_true$PFS_2_event <- ifelse( prog_true$Prog_II == 1 | prog_true$OS_event_death == 1, 1, 0)
prog_true$PFS_2_event %>% table
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
                                                      as.Date(prog_true$LAST_FULLOW_UP,origin="1970-01-01")))
# calulate the PFS_2_time
prog_true$PFS_2_time <- prog_true$PFS_2_date - prog_true$PFS_I_date
# calulate the PFS_2_months
prog_true$PFS_2_months <- round(prog_true$PFS_2_time / 30.5, 0) 
# FLAG
prog_true$flag <- ifelse(prog_true$PROTOCOLLO != "BO2005" & prog_true$PROTOCOLLO != "EMN02" & 
                           prog_true$TTP_I_months <= 24 & prog_true$Prog_I==1 &
                           prog_true$TP53_adj_D >= 1.9 & prog_true$MUT_p53_D_SUB_CLON ==0 , 
                         1, 0)

DF2 <- prog_true %>% filter(flag==0)

DF2$del_TP53_SNP_R <- ifelse(DF2$TP53_adj_R < 1.9, 1, 0)

#--- ADD PTS ID ----

DF1.1 <- left_join(DF1, pts_id )
DF1.1$pt_ID %>% unique() %>% length()
DF1 <- DF1.1 %>% filter(!is.na(pt_ID))


DF2.1 <- left_join(DF2, pts_id )
# DF2.1 %>% select(UPN, pt_ID) %>% View
DF2 <- DF2.1 %>% filter(!is.na(pt_ID))


#============== creation of databse for swimmer plot ===============

data_D <- DF1 %>% select(pt_ID, MUT_p53_D_SUB_CLON, del_TP53_SNP_D, PFS_I_months, PFS_I_event, OS_MESI, OS_event_death)

data_R <- DF2 %>% select(pt_ID, MUT_P53_R_SUB_CLON, del_TP53_SNP_R, PFS_2_event, PFS_2_months )

pts_all <- full_join(data_D, data_R, by = "pt_ID") %>% 
  select(pt_ID, MUT_p53_D_SUB_CLON, MUT_P53_R_SUB_CLON, 
         del_TP53_SNP_D, del_TP53_SNP_R, 
         PFS_I_months:OS_event_death, PFS_2_event,PFS_2_months)

 
# old  filter(UPN != "3581") #remove pts UPN 3581 because unvalutable SNP purity at Diagnosis (SNP218) missing data



pts_swimmer <- pts_all %>% filter(MUT_p53_D_SUB_CLON==1 | MUT_P53_R_SUB_CLON==1 | del_TP53_SNP_D ==1 | del_TP53_SNP_R==1) 

pts_swimmer$Second_PFS_months <- pts_swimmer$PFS_I_months + pts_swimmer$PFS_2_months

pts_swimmer.m = melt(pts_swimmer %>% select(pt_ID, OS_MESI, PFS_I_months, MUT_p53_D_SUB_CLON, MUT_P53_R_SUB_CLON, del_TP53_SNP_D, del_TP53_SNP_R),
             id.var=c("pt_ID","OS_MESI","PFS_I_months"))


pts_swimmer.m.D = melt(pts_swimmer %>% select(pt_ID, MUT_p53_D_SUB_CLON, del_TP53_SNP_D),
                       id.var=c("pt_ID"))
pts_swimmer.m.D$value[pts_swimmer.m.D$value==0] <- NA


pts_swimmer.m.R = melt(pts_swimmer %>% select(pt_ID, PFS_I_months, MUT_P53_R_SUB_CLON, del_TP53_SNP_R),
                       id.var=c("pt_ID", "PFS_I_months"))
pts_swimmer.m.R$valuetime <- ifelse(pts_swimmer.m.R$value == 1, pts_swimmer.m.R$PFS_I_months, NA)


# pts_swimmer.alive = pts_swimmer %>% select(UPN, OS_event_death, OS_MESI)




ORDER <- order(!is.na(pts_swimmer$MUT_P53_R_SUB_CLON) & !is.na(pts_swimmer$del_TP53_SNP_R), pts_swimmer$OS_MESI)

pts_swimmer$pt_ID <- factor(pts_swimmer$pt_ID, levels = pts_swimmer$pt_ID[ORDER], ordered = T)

pts_swimmer %>% ggplot(aes(pt_ID, OS_MESI)) +
  
  geom_bar(stat="identity", colour="black", fill="grey90", width=0.8)+
  
  # alterations D
  geom_point(data=pts_swimmer.m.D, 
             aes(pt_ID, -value, colour=variable, shape=variable), size=4) +
  
  # alterations R
  geom_point(data=pts_swimmer.m.R, 
             aes(pt_ID, valuetime, colour=variable, shape=variable), size=4) +
  
  # deaths
  geom_point(data=pts_swimmer %>% filter(OS_event_death==1),
             aes(pt_ID, OS_MESI+1), size=2.5, shape=4) +
  
  # alive
  geom_segment(data=pts_swimmer %>% filter(OS_event_death==0),
               aes(x=pt_ID, xend=pt_ID, y=OS_MESI + 0.5, yend=OS_MESI + 2),
               pch=25, size=0.5, arrow=arrow(type="closed", length=unit(0.07,"in"))) +
  
  # relapsed with EVAL
  geom_point(data=pts_swimmer %>% filter(PFS_I_event==1 & !is.na(MUT_P53_R_SUB_CLON) & !is.na(del_TP53_SNP_R)),
               aes(x=UPN, y= PFS_I_months-0.8), shape=124, size=3) +
  
  # relapsed no EVAL
  geom_point(data=pts_swimmer %>% filter(PFS_I_event==1 & (is.na(MUT_P53_R_SUB_CLON) | is.na(del_TP53_SNP_R))),
             aes(x=pt_ID, y= PFS_I_months-0.8), shape=1, size=3) +
  
  # SECOND relapsed no EVAL
  geom_point(data=pts_swimmer %>% filter(PFS_2_event==1),
             aes(x=pt_ID, y= Second_PFS_months-0.8), shape=1, size=3) +
  
  
  coord_flip() +
  
  scale_shape_manual(values = c(20,20,8,8) ) +
  # scale_fill_manual(values=hcl(seq(15,375,length.out=5)[1:4],100,70)) +
  # scale_colour_manual(values=c(hcl(seq(15,375,length.out=3)[1:2],100,40),"black")) +
  scale_y_continuous(limits=c(-1,172), breaks=seq(0,180,10), expand = expand_scale(mult = c(0.01,0))) +
  
  labs(y="Months") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())


######################### ONLY RELAPSED EVAL #######################

RE_pts_swimmer <- pts_swimmer %>% filter(!is.na(MUT_P53_R_SUB_CLON) & !is.na(del_TP53_SNP_R))
RE_pts_swimmer.m <- pts_swimmer.m %>% filter(pt_ID %in% unique(RE_pts_swimmer$pt_ID))
RE_pts_swimmer.m.D <- pts_swimmer.m.D %>% filter(pt_ID %in% unique(RE_pts_swimmer$pt_ID))
RE_pts_swimmer.m.R <- pts_swimmer.m.R %>% filter(pt_ID %in% unique(RE_pts_swimmer$pt_ID))


ORDER <- order( RE_pts_swimmer$OS_MESI)

RE_pts_swimmer$pt_ID <- factor(RE_pts_swimmer$pt_ID, levels = RE_pts_swimmer$pt_ID[ORDER], ordered = T)

RE_pts_swimmer %>% ggplot(aes(pt_ID, OS_MESI)) +
  
  geom_bar(stat="identity", colour="black", fill="grey90", width=0.8)+
  
  # relapsed with EVAL
  geom_point(data=RE_pts_swimmer %>% filter(PFS_I_event==1 & !is.na(MUT_P53_R_SUB_CLON) & !is.na(del_TP53_SNP_R)),
             aes(x=pt_ID, y= PFS_I_months-0.8), shape=05, size=3) +
  
  # relapsed no EVAL
  geom_point(data=RE_pts_swimmer %>% filter(PFS_I_event==1 & (is.na(MUT_P53_R_SUB_CLON) | is.na(del_TP53_SNP_R))),
             aes(x=pt_ID, y= PFS_I_months-0.8), shape=06, size=3) +
  
  # SECOND relapsed no EVAL
  geom_point(data=RE_pts_swimmer %>% filter(PFS_2_event==1),
             aes(x=pt_ID, y= Second_PFS_months-0.8), shape=06, size=3) +
  
  # alterations D
  geom_point(data=RE_pts_swimmer.m.D, 
             aes(pt_ID, -value, colour=variable, shape=variable), size=4) +
  
  # alterations R
  geom_point(data=RE_pts_swimmer.m.R, 
             aes(pt_ID, valuetime, colour=variable, shape=variable), size=4) +
  
  # deaths
  geom_point(data=RE_pts_swimmer %>% filter(OS_event_death==1),
             aes(pt_ID, OS_MESI+2), size=2.5, shape=4) +
  
  # alive
  geom_segment(data=RE_pts_swimmer %>% filter(OS_event_death==0),
               aes(x=pt_ID, xend=pt_ID, y=OS_MESI + 0.5, yend=OS_MESI + 3),
               pch=25, size=0.5, arrow=arrow(type="closed", length=unit(0.07,"in"))) +
  
  coord_flip() +
  
  scale_colour_discrete(name = "Event", labels = c("del TP53 Diagnosis", "del TP53 Relapse", "mut TP53 Diagnosis","mut TP53 Relapse")) +
  scale_shape_manual(values = c(20,20,8,8), name = "Event", labels = c("del TP53 Diagnosis", "del TP53 Relapse", "mut TP53 Diagnosis","mut TP53 Relapse") ) +
  
  # scale_fill_manual(values=hcl(seq(15,375,length.out=5)[1:4],100,70)) +
  # scale_colour_manual(values=c(hcl(seq(15,375,length.out=3)[1:2],100,40),"black")) +
  scale_y_continuous(limits=c(-1,172), breaks=seq(0,180,10), expand = expand_scale(mult = c(0.01,0))) +
  
  labs(y="Months") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())


dummy <- data.frame(x=c(1,2,3,4),y=c(1,2,3,4),event=c("death","alive","relapse analyzed", "relapse not analyzed"))

dummy %>% ggplot(aes(x,y)) + geom_point(aes(shape=event)) + 
  scale_shape_manual(values = c(4,95,5,6), name="Event", labels=c("death","alive","relapse analyzed", "relapse not analyzed")) +
  theme_bw() 
