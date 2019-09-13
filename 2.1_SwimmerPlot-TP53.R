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
