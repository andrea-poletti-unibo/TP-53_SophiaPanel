```{r}
library(tidyverse)
library(RODBC)

db <- odbcConnectAccess2007("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT SophiaPanel - TP53 - Documenti/TP53_DB_v1.accdb")

df <- sqlFetch(db, "Query_Survival_Analysis")
```

```{r}
df$PROTOCOLLO %>% table
```

```{r}
amb <- df %>% filter(PROTOCOLLO != "BO2005" & PROTOCOLLO != "EMN02")

amb$TERAPIA_INDUZIONE_I_LINEA %>% table
```


```{r}

amb24 <- amb %>% filter(TTP_I_months <= 24 & Prog_I==1 )

amb24$TERAPIA_INDUZIONE_I_LINEA %>% as.character() %>% table

```


```{r}
table(amb24$TP53_adj < 1.9, amb24$MUT_p53_D_SUB_CLON)


df$flag <- ifelse(df$PROTOCOLLO != "BO2005" & df$PROTOCOLLO != "EMN02" & 
                    df$TTP_I_months <= 24 & df$Prog_I==1 &
                    df$TP53_adj >= 1.9 & df$MUT_p53_D_SUB_CLON ==0 , 
                  1, 0)

```



```{r include=FALSE}
DF <- df %>% filter(flag != 1)


DF$call10_del_TP53 <- ifelse(DF$TP53_adj <= 1.9, 1,0)
DF$call10_broad_amp_1q <- ifelse(DF$broad_1q >= 2.1, 1, 0)
DF$call10_broad_del_1p <- ifelse(DF$broad_1p <= 1.9, 1, 0)
DF$call10_broad_del_13q <- ifelse(DF$broad_13q <= 1.9, 1, 0)
DF$call10_broad_del_17p <- ifelse(DF$broad_17p <= 1.9, 1, 0)
```


#============ Survival analysis =====================
```{r include=FALSE}
library(data.table)
library(tidyverse)
library(survival)
library(survminer)


OS <- Surv( time = DF$OS_MESI, event = DF$OS_event_death)
PFS <- Surv( time = DF$PFS_I_months, event = DF$PFS_I_event)
```

# TP53 -> significative in OS and PFS
```{r}
ggsurvplot(survfit(OS ~ DF$call10_del_TP53, data = DF), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ DF$call10_del_TP53, data = DF), pval = T, risk.table = T, xlab = "PFS")


coxph(OS ~ DF$call10_del_TP53 + strata(DF$ISS), data = DF) %>% summary
coxph(PFS ~ DF$call10_del_TP53 + strata(DF$ISS) , data = DF) %>% summary


```

# 1q -> only OS
```{r}

sum(DF$call10_broad_amp_1q)
61/143

ggsurvplot(survfit(OS ~ DF$call10_broad_amp_1q, data = DF), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ DF$call10_broad_amp_1q, data = DF), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ DF$call10_broad_amp_1q + strata(ISS), data = DF) %>% summary
coxph(PFS ~ DF$call10_broad_amp_1q + strata(ISS), data = DF) %>% summary

```

# 1p -> only OS
```{r}
sum(DF$call10_broad_del_1p)
35/143

ggsurvplot(survfit(OS ~ DF$call10_broad_del_1p, data = DF), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ DF$call10_broad_del_1p, data = DF), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ DF$call10_broad_del_1p + strata(ISS), data = DF) %>% summary
coxph(PFS ~ DF$call10_broad_del_1p + strata(ISS), data = DF) %>% summary
```

# 13q -> only OS
```{r}
sum(DF$call10_broad_del_13q)
78/143

ggsurvplot(survfit(OS ~ DF$call10_broad_del_13q, data = DF), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ DF$call10_broad_del_13q, data = DF), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ DF$call10_broad_del_13q + strata(ISS), data = DF) %>% summary
coxph(PFS ~ DF$call10_broad_del_13q + strata(ISS), data = DF) %>% summary
```

# t 4 14 -> no sig
```{r}
sum(DF$Fish_T_4_14, na.rm = T)
29/143

ggsurvplot(survfit(OS ~ DF$Fish_T_4_14, data = DF), pval = T, risk.table = T, xlab = "OS")
ggsurvplot(survfit(PFS ~ DF$Fish_T_4_14, data = DF), pval = T, risk.table = T, xlab = "PFS")

coxph(OS ~ DF$Fish_T_4_14 + strata(ISS), data = DF) %>% summary
coxph(PFS ~ DF$Fish_T_4_14 + strata(ISS), data = DF) %>% summary
```


#======= multivariate =========

# OS: TP53 + amp 1q + del 13 + del 1p
```{r}
coxph(OS ~ DF$call10_del_TP53 + DF$call10_broad_amp_1q + DF$call10_broad_del_13q + DF$call10_broad_del_1p + strata(ISS), data = DF) %>% summary
```

# PFS: TP53 + del 13
```{r}
coxph(PFS ~ DF$call10_del_TP53 + DF$call10_broad_del_13q + strata(ISS), data = DF) %>% summary

```

