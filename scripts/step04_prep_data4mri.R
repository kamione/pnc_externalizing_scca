# Information ------------------------------------------------------------------
# Version    : 0.2.0
# Date       : 30-March-2020
# Description: Factor analysis and calculate the factor scores 

# Working Directory ------------------------------------------------------------
path_proj <- "/Users/one_macbook/Desktop"
path_wd   <- file.path(path_proj, "2020_pnc_fc_extbeh")
setwd(path_wd)

# Load Packages ----------------------------------------------------------------
library(tidyverse)
source(file.path("utilities", "R", "customized_functions.R"))

# Data I/O ---------------------------------------------------------------------
path_mriinfo <- file.path("data", "raw", "PNC_MRIQuality.xlsx")

mri_info <- path_mriinfo %>%
  read_excel(sheet = 1) %>%
  rename(SUBJID = SubjID) %>%
  filter(T1 == "Usable") %>%
  filter(`no. of img` == 124) %>%
  filter(`Mean FD` < 0.2) %>%
  rename(MeanFD = `Mean FD`) %>%
  select(SUBJID, MeanFD)

# visualize whole samples (n = 8,832) ADHD ~ DBD
ext_sum_whole <- read_rds(file.path("outputs", "psychopathology_sum.rds")) %>%
  select(SUBJID, ADD, CDD, ODD) %>%
  mutate(DBD = CDD + ODD) %>%
  select(SUBJID, ADD, DBD) %>%
  mutate(ADD = scale_range01(ADD)) %>%
  mutate(DBD = scale_range01(DBD))

rval_whole <-cor(ext_sum_whole$ADD, ext_sum_whole$DBD)
p1_hexbin <- ext_sum_whole %>%
  ggplot(aes(x = ADD, y = DBD)) +
  geom_hex(bins = 9) +
  stat_smooth(color = "tomato3", method = 'glm') +
  scale_fill_gradientn(colours = c("darkslateblue", "cyan3"), trans = "log10") +
  guides(fill = guide_legend(title = "Frequency", title.position = "left",
                             direction = "horizontal")) + 
  annotate("text", x = 0, y = 1.05, color = "gray10", size = 5, hjust = 0.2,
           label = sprintf("italic(r) == %.3f", rval_whole), parse = TRUE) +
  labs(x = "ADHD Symptom Loads", y = "DBD Symptom Loads") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.1, 0)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0, 0.1)) +
  theme_classic() +
  theme(legend.position = "bottom")

# visualize MRI samples (n = 1,149) ADHD ~ DBD
ext_sum_mri <- mri_info %>%
  select(SUBJID) %>%
  merge(ext_sum_whole, by = "SUBJID")

rval_mri <-cor(ext_sum_mri$ADD, ext_sum_mri$DBD)
p2_hexbin <- ext_sum_mri %>%
  ggplot(aes(x = ADD, y = DBD)) +
  geom_hex(bins = 9) +
  stat_smooth(color = "tomato3", method = 'glm') +
  scale_fill_gradientn(colours = c("darkslateblue", "cyan3"), trans = "log10") +
  guides(fill = guide_legend(title = "Frequency", title.position = "left",
                             direction = "horizontal")) + 
  annotate("text", x = 0, y = 1.05, color = "gray10", size = 5, hjust = 0.2,
           label = sprintf("italic(r) == %.3f", rval_mri), parse = TRUE) +
  labs(x = "ADHD Symptom Loads", y = "DBD Symptom Loads") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.1, 0)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0, 0.1)) +
  theme_classic() +
  theme(legend.position = "bottom")

hexbin_plot <- ggarrange(p1_hexbin, p2_hexbin, nrow = 1, ncol = 2)
ggexport(hexbin_plot, filename = file.path("figures", "ext_hexbinplot.pdf"), width = 9, height = 4.5)

# write the preprocessed data to csv file
mri_data <- read_rds(file.path("outputs", "imputated_data.rds")) %>%
  select(SUBJID:Race3, WRAT_CR_RAW) %>%
  rename(Age = ageAtClinicalAssess) %>%
  rename(`WRAT IQ` = WRAT_CR_RAW) %>%
  merge(mri_info, by = "SUBJID") %>%
  merge(ext_sum_mri, by = "SUBJID")

write_csv(mri_data, file.path("data", "preprocessed", "pnc_n1149_data"))




