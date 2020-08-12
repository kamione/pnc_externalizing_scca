# Information ------------------------------------------------------------------
# Version    : 0.2.0
# Date       : 30-March-2020
# Description: Imputating the missing scores in symptoms and CNB based on age,
#              sex and race

# Working Directory ------------------------------------------------------------
# please change this folder directory
path_proj <- "/Users/one_macbook/Desktop"
path_wd   <- file.path(path_proj, "2020_pnc_fc_extbeh")
setwd(path_wd)

# Load Packages ----------------------------------------------------------------
library(tidyverse)
library(readxl)
library(ggplot2)
library(ggpubr)
library(mice)

# Data I/O ---------------------------------------------------------------------
# path to file directories
path_psy_items <- file.path("data", "raw", "psychoapthology_nonskip_items.csv")
path_cnb_items <- file.path("data", "raw", "CNB_items.csv")
path_phenotype <- file.path("data", "raw", "phs000607.v3.pht003445.v3.p2.c1.Neurodevelopmental_Genomics_Subject_Phenotypes.GRU-NPU.xlsx")
path_demogr    <- file.path("data", "raw", "phs000607.v3.pht007881.v1.p2.c1.Interview_Dates_Update.GRU-NPU.txt")

# CNB Abbreivations
# PCET: Penn Conditional Exclusion Test
# PCPT: Penn Continuous Performance Test
# LNB : Letter N-Back Task
# PWMT: Penn Word Memory Task
# PFMT: Penn Face Memory Task
# VLOT: Visual Object Learning Test
# PVRT: Penn Verbal Reasoning Test
# PMRT: Penn Matrix Reasoning Test
# PLOT: Penn Line Orientation Test
# PEIT: Penn Emotion Identification Test
# PEDT: Penn Emotion Differentiation Test
# PADT: Penn Age Differentiation Test
# WART: Wide Range Assessment Test (Reading IQ)

psy_items <- read_csv(path_psy_items) %>% pull(Item)
cnb_items <- read_csv(path_cnb_items) %>% pull(Item)

phenotype <- path_phenotype %>% 
  read_excel(sheet = 1, na = c("NA", " ", 9)) %>%
  select(c(SUBJID, INT_TYPE, Sex, Race, psy_items, cnb_items))
age <- path_demogr %>% read_tsv(skip = 10) %>%
  select(c(SUBJID, ageAtClinicalAssess))

# Data Merging -----------------------------------------------------------------
# 9,337 participants
data <- merge(age, phenotype, by = "SUBJID") %>%
  transform(SUBJID = as.character(SUBJID)) %>%
  subset(INT_TYPE == "AP" | INT_TYPE == "MP" | INT_TYPE == "YPI") %>%
  select(-INT_TYPE) %>%
  filter(Sex == "M" | Sex == "F") %>%
  subset(Race != "NA") %>%
  mutate(Sex = as.factor(recode(Sex, "M" = 1, "F" = 0))) %>%
  mutate(Race = recode(Race, "EA" = 1, "AA" = 2, .default = 3)) %>%
  mutate(one = 1) %>%
  spread(Race, one, fill = 0, sep = "") %>%
  select(SUBJID:Sex, Race1:Race3, psy_items, cnb_items)

# removed subjects with more than 10 missing variables based on histogram
hist_fig <- data %>%
  is.na() %>%
  rowSums() %>%
  hist(col = "gray30", xlab = "Numbers of Missings", main = "")

data_filtered <- data[data %>% is.na() %>% rowSums() <= 10, ]

# visualize missing datq
missing_values <- data_filtered %>%
  gather(key = "key", value = "val") %>%
  mutate(isna = is.na(val)) %>%
  group_by(key) %>%
  mutate(total = n()) %>%
  group_by(key, total, isna) %>%
  summarise(num.isna = n()) %>%
  mutate(pct = num.isna / total) 

levels <- (missing_values %>% filter(isna == TRUE) %>% arrange(desc(pct)))$key

proportion_plot <- missing_values %>%
  ggplot() +
  geom_bar(aes(x = reorder(key, desc(pct)), y = pct, fill=isna), stat = 'identity', alpha = 0.8) +
  scale_x_discrete(limits = levels) +
  scale_fill_manual(name = "", values = c('grey20', 'tomato2'), labels = c("Present", "Missing")) +
  coord_flip() +
  theme(axis.text.y = element_text(size = 4)) +
  theme(legend.position = "none") +
  labs(x = "Variable", y = "Proportion of Missings") 
 
row_plot <- data_filtered %>%
  mutate(id = row_number()) %>%
  gather(-id, key = "key", value = "val") %>%
  mutate(isna = is.na(val)) %>%
  ggplot(aes(key, id, fill = isna)) +
  geom_raster(alpha = 0.8) +
  scale_fill_manual(name = "",
                    values = c('grey20', 'tomato2'),
                    labels = c("Present", "Missing")) +
  scale_x_discrete(limits = levels) +
  labs(y = "Missings in Individuals") +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_text(size = 4)) +
  theme(legend.position = "none") +
  coord_flip()

fig_missingplot <- ggarrange(proportion_plot, row_plot, ncol = 2, nrow = 1)
ggexport(fig_missingplot, filename = "figures/missingplot.pdf", width = 10, height =  12)

# Imputatoin -------------------------------------------------------------------
# using "mice" packages, imputated missings based on age, sex and race
imputated_scores <- data_filtered %>% select(psy_items, cnb_items)
n_items <- dim(imputated_scores)[2]
for (ith_item in 1:n_items) {
  item <- colnames(imputated_scores[ith_item])
  input <- data_filtered %>% select(ageAtClinicalAssess:Race3, item)
  n_var <- dim(input)[2]
  if (input[n_var] %>% is.na() %>% sum() == 0) {next}
  # flag binary target
  response_set <- input[, n_var] %>% unique() %>% discard(is.na) %>% sort()
  # imputate based on condition
  if (identical(response_set, c(0, 1))) {
    imputation <- mice(input, m = 5, maxit = 50, defaultMethod = "logreg",
                       seed = 1234, printFlag = FALSE)
    imputated_scores[ith_item] <- complete(imputation, 1)[n_var]
  } else if (!identical(response_set, c(0, 1)) && length(response_set) > 1) {
    imputation <- mice(input, m = 5, maxit = 50, defaultMethod = "pmm",
                       seed = 1234, printFlag = FALSE)
    imputated_scores[ith_item] <- complete(imputation, 1)[n_var]    
  } else if (length(response_set) == 1) {
    imputated_scores[ith_item] <- response_set
  }
}

# replace missing with imputated scores
imputated_data <- data_filtered
imputated_data[c(psy_items, cnb_items)] <- imputated_scores

# double check
flag <- imputated_data %>% is.na() %>% any()
if (flag == 1) {statement = "Yes"} else {statement = "No"}
msg  <- paste("Any missings:", statement, sep = " ")
print(msg)

# save imputated dataframe as .rds
saveRDS(imputated_data, "outputs/imputated_data.rds")
