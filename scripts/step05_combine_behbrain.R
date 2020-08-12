# Information ------------------------------------------------------------------
# Version    : 0.1.0
# Date       : 13-April-2020
# Description: prepare data for sCCA

# Working Directory ------------------------------------------------------------
# please change to the working directory
path_wd <- "/Users/one_macbook/Desktop/2020_pnc_fc_extbeh"
setwd(path_wd)

# Load Library -----------------------------------------------------------------
library(tidyverse)
library(caret)

# Data I/O ---------------------------------------------------------------------
# indicate paths
path_behdata   <- file.path('data', 'preprocessed', 'pnc_n945_data.csv')    
path_fconn_wb  <- file.path('outputs', 'shen_fconn_selected.csv')
path_fconn_mcc <- file.path('outputs', 'mcc_fcmaps_shen.csv')
path_fconn_amy <- file.path("outputs", "amy_fcmaps_shen.csv")
path_fconn_striatum <- file.path("outputs", "striatum_fcmaps_shen.csv")
path_fconn_ipl <- file.path("outputs", "ipl_fcmaps_shen.csv")
path_fconn_presma <- file.path("outputs", "presma_fcmaps_shen.csv")

#path_fconn_mcc <- file.path('outputs', 'mcc_fcmaps_selected.csv')
path_subjlist1 <- file.path('outputs', 'wholebrain_subjlist.txt')
path_subjlist2 <- file.path('outputs', 'fcmaps_subjlist.txt')

# read in data
fconn_wb_data <- read_tsv(path_subjlist1, col_names = "SUBJID") %>%
  cbind(read_csv(path_fconn_wb, col_names = FALSE))
fconn_mcc_data <- read_tsv(path_subjlist2, col_names = "SUBJID") %>% 
  cbind(read_csv(path_fconn_mcc, col_names = FALSE))
fconn_amy_data <- read_tsv(path_subjlist2, col_names = "SUBJID") %>% 
  cbind(read_csv(path_fconn_amy, col_names = FALSE))
fconn_striatum_data <- read_tsv(path_subjlist2, col_names = "SUBJID") %>% 
  cbind(read_csv(path_fconn_striatum, col_names = FALSE))
fconn_ipl_data <- read_tsv(path_subjlist2, col_names = "SUBJID") %>% 
  cbind(read_csv(path_fconn_ipl, col_names = FALSE))
fconn_presma_data <- read_tsv(path_subjlist2, col_names = "SUBJID") %>% 
  cbind(read_csv(path_fconn_presma, col_names = FALSE))
#fconn_sma <- read_csv(path_fconn_sma, col_names = FALSE)
#fconn_mcc <- read_csv(path_fconn_mcc, col_names = FALSE)
#fconn_mcc <- read_csv(path_fconn_mcc, col_names = FALSE)
#fconn_mcc <- read_csv(path_fconn_mcc, col_names = FALSE)

data_wb <- path_behdata %>%
  read_csv() %>%
  merge(fconn_wb_data)

oldnames <- paste("X", 1:268, sep = "")

data_mcc <- path_behdata %>%
  read_csv() %>%
  merge(fconn_mcc_data) %>% 
  rename_at(vars(oldnames), ~ paste("mcc", 1:268, sep = ""))

data_amy <- path_behdata %>%
  read_csv() %>%
  merge(fconn_amy_data) %>% 
  rename_at(vars(oldnames), ~ paste("amy", 1:268, sep = ""))

data_striatum <- path_behdata %>%
  read_csv() %>%
  merge(fconn_striatum_data) %>% 
  rename_at(vars(oldnames), ~ paste("striatum", 1:268, sep = ""))

data_ipl <- path_behdata %>%
  read_csv() %>%
  merge(fconn_ipl_data) %>% 
  rename_at(vars(oldnames), ~ paste("ipl", 1:268, sep = ""))

data_presma <- path_behdata %>%
  read_csv() %>%
  merge(fconn_presma_data) %>% 
  rename_at(vars(oldnames), ~ paste("presma", 1:268, sep = ""))


# split the data into discovery and replication samples
set.seed(1234)
tridx <- createDataPartition(data$Age, p = 0.667, list = FALSE, times = 1)

data_wb_discovery   <- data_wb[tridx, ]
data_wb_replication <- data_wb[-tridx, ]
saveRDS(data_wb, file.path("data", "preprocessed", "data_wb.rds"))
saveRDS(data_wb_discovery, file.path("data", "preprocessed", "data_wb_discovery.rds"))
saveRDS(data_wb_replication, file.path("data", "preprocessed", "data_wb_replication.rds"))

data_mcc_discovery   <- data_mcc[tridx, ]
data_mcc_replication <- data_mcc[-tridx, ]
saveRDS(data_mcc, file.path("data", "preprocessed", "data_mcc_shen.rds"))
saveRDS(data_mcc_discovery, file.path("data", "preprocessed", "data_mcc_discovery.rds"))
saveRDS(data_mcc_replication, file.path("data", "preprocessed", "data_mcc_replication.rds"))


saveRDS(data_wb, file.path("data", "preprocessed", "data_wb_shen.rds"))
saveRDS(data_mcc, file.path("data", "preprocessed", "data_mcc_shen.rds"))
saveRDS(data_amy, file.path("data", "preprocessed", "data_amy_shen.rds"))
saveRDS(data_striatum, file.path("data", "preprocessed", "data_striatum_shen.rds"))
saveRDS(data_ipl, file.path("data", "preprocessed", "data_ipl_shen.rds"))
saveRDS(data_presma, file.path("data", "preprocessed", "data_presma_shen.rds"))
