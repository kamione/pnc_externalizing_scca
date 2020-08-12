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
library(mirt)
library(psych)
library(pheatmap)
library(ggplot2)
library(wesanderson)
library(readxl)
library(corrplot)

source(file.path("utilities", "R", "customized_functions.R"))

# Load Imputated Data ----------------------------------------------------------
imputated_data <- readRDS(file.path("outputs", "imputated_data.rds"))
path_age       <- "data/phs000607.v3.pht007881.v1.p2.c1.Interview_Dates_Update.GRU-NPU.txt"
path_cnb_items <- "data/CNB_items.csv"

# Psychopathology Factor Analysis ----------------------------------------------
psy_items <- imputated_data %>% select(ADD011:SUI002)

psychopathology_sum <- psy_items %>%
  colnames() %>%
  str_extract("...") %>%
  unique() %>%
  sapply(function(x) rowSums(psy_items[, grep(x, names(psy_items)), drop = FALSE]))

n_labels <- dim(psychopathology_sum)[2]
psy_labels <- colnames(psychopathology_sum)
zerodistribution <- data.frame(
  "label"  = character(n_labels),
  "n_zero" = double(n_labels),
  stringsAsFactors = FALSE
)

for (ith_label in 1:n_labels) {
  tmp <- table(psychopathology_sum[, ith_label])
  if (is.null(tmp[["0"]]) == 1) {
    zerodistribution[ith_label, ]  = c(psy_labels[ith_label], 0)
  } else {
    zerodistribution[ith_label, ]  = c(psy_labels[ith_label], tmp[["0"]])
  }
}
zerodistribution$n_zero <- as.numeric(zerodistribution$n_zero)

# visualize no. of zero cases in each diagnostic label
zerodistribution_plot <- zerodistribution %>%
  ggplot(aes(x = .[, 1], y = .[, 2])) + 
    geom_bar(stat="identity") +
    labs(x = "Common Psychopathology",
         y = "Number of People Showing Null Symptom") +
    scale_y_continuous(breaks = seq(0, 7000, 1000), expand = c(0, 100)) +
    theme_classic()
ggsave(file.path("figures", "zero_distribution.pdf"), zerodistribution_plot, width = 6, height = 4)

# visualize common mental health disorder correlation in PNC
psy_corr    <- cor(psychopathology_sum)

color_break <- c("gray40", "white", "cyan4", "yellow", "indianred2")
color_map   <- colorRampPalette(color_break)(256)
pdf(file.path("figures", "psychoapthology_correlation.pdf"), width = 8, height = 8)
fig_psycorr <- psy_corr %>%
  corrplot(method = "color",  col = color_map,
           type = "lower",
           addCoef.col = "gray30", # Add coefficient of correlation
           tl.col = "gray30", tl.srt = 30, #Text label color and rotation
           diag = FALSE,
           cl.lim = c(0, 1)
  )
dev.off()

imputated_data_with_psysum <- bind_cols(imputated_data[1], data.frame(psychopathology_sum))
saveRDS(imputated_data_with_psysum, file.path("outputs", "psychopathology_sum.rds"))

# determine number of factors based on parallel analysis
set.seed(1234)
parallel <- fa.parallel(psychopathology_sum, cor = "cor", fa = "fa",
                        m = "minres", n.iter = 1000, quant = .95, plot = FALSE)

screeplot(parallel) # factor 5, 6 and 7 are similar
ggsave(file.path("figures", "pa_screeplot.png"), width = 5, height = 4)

# run factor 5 to 7
set.seed(1234)
mod5 <- mirt(psychopathology_sum, model = 5, itemtype = "graded", method = "MHRM")
bf5_loadings <- summary(mod5, rotate = "bifactorT", suppress = 0.1)
set.seed(1234)
mod6 <- mirt(psychopathology_sum, model = 6, itemtype = "graded", method = "MHRM")
bf6_loadings <- summary(mod6, rotate = "bifactorT")
set.seed(1234)
mod7 <- mirt(psychopathology_sum, model = 7, itemtype = "graded", method = "MHRM")
bf7_loadings <- summary(mod7, rotate = "bifactorT", suppress = 0.1)

# compare the models
diff1 <- anova(mod5, mod6)
diff2 <- anova(mod5, mod7)
##  model 5 (1 + 4) is a better model consistent to the literature

# adjust bifacotr loadings (some are reversed)
bf5_loadings_adjusted <- bf5_loadings
bf5_loadings_adjusted$rot[, c(1, 5)] <- -bf5_loadings$rot[, c(1, 5)]

psy_fa_heatmap <- bf5_loadings_adjusted$rot %>%
  pheatmap(
    fontsize        = 12,
    treeheight_row  = 0,
    cluster_cols    = FALSE,
    legend_breaks   = seq(-0.2, 0.6, by = 0.2),
    display_numbers = TRUE,
    fontsize_number = 10,
    labels_col      = c("Psychopathology", "Fear", "Externalizing", "Anxious/Misery", "Psychosis"),
    angle_col       = 90,
    filename        = file.path("figures", "psych_fa_heatmap.pdf"),
    width           = 5,
    height          = 9,
    silent          = TRUE
)

fscores_psy_5fa <- fscores(mod5, rotate = "bifactorT", QMC = TRUE)
fscores_psy_5fa_adjusted <- fscores_psy_5fa
fscores_psy_5fa_adjusted[, c(1, 4, 5)] <- -fscores_psy_5fa[, c(1, 4, 5)]
psy_fscores_data <- bind_cols(imputated_data[1], data.frame(fscores_psy_5fa_adjusted))
colnames(psy_fscores_data) <- c("SUBJID", "Psychopathology", "Fear", "Externalizing", "Anxious/Misery", "Psychosis")
saveRDS(psy_fscores_data, file.path("outputs", "psy_fscores_data.rds"))

# CNB Factor Analysis ----------------------------------------------------------
cnb_items <- imputated_data %>% select(PCET_ACC2:WRAT_CR_RAW)
task_names <- read_csv(path_cnb_items) %>% pull(Task_Abbreviations)
task_names4var <- paste(task_names[1:12], "efficiency", sep = "_")

cnb_items_zscore <- as.data.frame(scale(cnb_items))
for (ith_var in 1:12) {
  cnb_items_zscore[[task_names4var[ith_var]]] <- cnb_items_zscore %>%
    select(starts_with(task_names[ith_var])) %>%
    reduce(`-`)
}

cnb_efficiency <- cnb_items_zscore %>% select(task_names4var)

# run unidimenstional factor analysis
set.seed(1234)
mod_cnb <-  fa(cnb_efficiency)
# compare the models
fscores_cnb_1fa <- mod_cnb$scores

fscores_cnb_data <- bind_cols(imputated_data[1], data.frame(fscores_cnb_1fa))
colnames(fscores_cnb_data)[2] <- "General Cognition"

# visualize correlation of all facotr items
all <- cbind(fscores_psy_5fa_adjusted, fscores_cnb_1fa, cnb_items[, 25])
colnames(all) <- c("Psychopathology", "Fear", "Externalizing", "Anxious/Misery",
                   "Psychosis", "General Cognition", "WRAT IQ")
all_corr <- cor(all)

color_break <- c("steelblue", "cyan4", "white", "yellow", "indianred2")
color_map   <- colorRampPalette(color_break)(256)
pdf(file.path("figures", "factor_correlation.pdf"), width = 8, height = 8)
fig_facorr <- all_corr %>%
  corrplot(method = "color",  col = color_map,
           type = "lower",
           addCoef.col = "gray30", # Add coefficient of correlation
           tl.col = "gray30", tl.srt = 45, #Text label color and rotation
           diag = FALSE
  )
dev.off()

saveRDS(fscores_cnb_data, file.path("outputs", "fscores_cnb_data.rds"))




