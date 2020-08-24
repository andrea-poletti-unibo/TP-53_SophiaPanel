# library(RODBC) 
library(odbc)
library(DBI)
library(data.table)
library(tidyverse)

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

df2 <- df %>% filter(flag != 1)

df2$pt_ID <- paste0("ID-", sprintf("%03d", 1:143) ) 

export <- df2 %>% select(UPN, pt_ID)

write_tsv(export, "valid_pts_ID_table.txt")

