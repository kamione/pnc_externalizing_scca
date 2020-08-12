# Information ------------------------------------------------------------------
# Version    : 0.1.0
# Date       : 13-April-2020
# Description: sCCA to link connectome and behavior 

# Working Directory ------------------------------------------------------------
# please change to the working directory
path_wd <- "/Users/one_macbook/Desktop/2020_pnc_fc_extbeh"
setwd(path_wd)

# Load Packages ----------------------------------------------------------------
library(tidyverse)
library(readxl)
library(parallel)
library(caret)
library(PMA)
library(matrixStats)
library(psych)
library(plyr)
library(gghighlight)
library(scales) # adjust color alpha
library(lattice) # leveloplot
library(viridis)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(fmsb) # radarchart
library(scales) # scale alpha value of color
source(file.path("utilities", "R", "scca_functions.R"))
source(file.path("utilities", "R", "customized_scca_functions.R"))


# Data I/O ---------------------------------------------------------------------
# get data
data_mcc       <- readRDS(file.path("data", "preprocessed", "data_mcc_shen.rds"))
fconn_mcc      <- data_mcc %>% select(starts_with("mcc"))
data_presma    <- readRDS(file.path("data", "preprocessed", "data_presma_shen.rds"))
fconn_presma   <- data_presma %>% select(starts_with("presma"))
data_striatum  <- readRDS(file.path("data", "preprocessed", "data_striatum_shen.rds"))
fconn_striatum <- data_striatum %>% select(starts_with("striatum"))

fconn <- cbind(fconn_mcc, fconn_presma, fconn_striatum)
ext   <- data_mcc %>% select(ADD011:ODD006)
data  <- cbind(data_mcc %>% select(SUBJID:DBD), fconn)

ext <- ext[-128, ]
fconn <- fconn[-128, ]
data <- data[-128, ]
data <- data %>% mutate(EXT = ADD + DBD)

# Data Preprocessing -----------------------------------------------------------
# regress out covariates
fconn_resid <- apply(
  fconn, 2, function(x) {
    glm(x ~ Age + as.factor(Sex) + as.factor(Race) + MeanFD, data = data) %>%
      residuals.glm(type = "response")
  }
)
ext_resid <- apply(
  ext, 2, function(x) {
    glm(x ~ Age + as.factor(Sex) + as.factor(Race),
        family = binomial(link = 'logit'), data = data) %>%
      residuals.glm(type = "response")
  }  
)

# Data Partition ---------------------------------------------------------------
set.seed(1111)
tridx <- createDataPartition(data$Age, p = 0.667, list = FALSE)
data_train <- NULL
data_train$fconn <- fconn_resid[tridx, ]
data_train$ext  <- ext_resid[tridx, ]
data_test <- NULL
data_test$fconn <- fconn_resid[-tridx, ]
data_test$ext  <- ext_resid[-tridx, ]

data_tr <- data[tridx, ]
data_te <- data[-tridx, ]


# sCCA on Discovery ------------------------------------------------------------
# create 3 folds 100 times in the sub-training set to obtain average
set.seed(1111)
resample_idx_tr     <- createDataPartition(data_tr$EXT, p = 0.667, list = TRUE, times = 100)
fconn_10resample_tr <- mclapply(resample_idx_tr, function(x) data_train$fconn[x,])
ext_10resample_tr   <- mclapply(resample_idx_tr, function(x) data_train$ext[x,])

pen_x <- seq(0.1, 1, 0.1)
pen_z <- 0.5
set.seed(1111)
scca_gs <- ccaDWfoldgs(fconn_10resample_tr, ext_10resample_tr, pen_x, pen_z)
scca_gs
pen_x_selected <- scca_gs$PENX
pen_z_selected <- 0.5

nItems <- 22

scca_main_tr <- run_scca(data_train$fconn, data_train$ext, pen_x_selected, pen_z_selected)

EVs <- scca_main_tr$d ^ 2 / sum(scca_main_tr$d ^ 2)
EVs_reorder <- sort(EVs, decreasing = TRUE, index.return = TRUE)
EVs_df <- data.frame(Mode = as.factor(1:nItems), EV = EVs_reorder$x)
expectedEV <- 1 / nItems
maxEV <- ceiling(max(EVs) * 100) / 100 + 0.005

EVs_screeplot <- EVs_df %>%
  mutate(highlight_flag = ifelse(Mode %in% 1:2, TRUE, FALSE)) %>% 
  ggplot(aes(Mode, EV, color = highlight_flag)) +
    geom_point(stat = 'identity', aes(size = EV)) +
    scale_color_manual(values = c("grey80", "grey35")) +
    #geom_hline(yintercept = expectedEV, linetype = "dashed", color = "grey50") +
    scale_x_discrete(breaks = seq(1, 22, 1)) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, maxEV),
      breaks = seq(0, maxEV, 0.02)
    ) +
    labs(x = "Mode", y = "Explained Variance") +
    theme_classic(base_size = 12) +
    theme(legend.position = 'none')
EVs_screeplot
ggsave(file.path("figures", "adhd_scca_train_screeplot.png"), EVs_screeplot, width = 6, height = 4, dpi = 400)

# run sCCA with permutation
nMode <- 3
scca_cand_tr <- run_scca(data_train$fconn, data_train$ext, pen_x_selected, pen_z_selected, nMode)

# create 3 folds 100 times in the sub-training set to obtain average
set.seed(1111)
resample_idx_tr     <- createDataPartition(data_tr$EXT, p = 0.667, list = TRUE, times = 100)
fconn_10resample_tr <- mclapply(resample_idx_tr, function(x) data_train$fconn[x,])
ext_10resample_tr   <- mclapply(resample_idx_tr, function(x) data_train$ext[x,])

scca_cca_tr  <- mclapply(
  seq_along(resample_idx_tr), function(i) {
    run_scca(fconn_10resample_tr[[i]], ext_10resample_tr[[i]], pen_x_selected, pen_z_selected, nMode)
  }
)
scca_cca_tr_ro     <- sapply(scca_cca_tr, function(x) reorderCCA(x, scca_cand_tr, nMode))
scca_cca_tr_cor    <- rowMeans(simplify2array(scca_cca_tr_ro['cors', ]), na.rm = TRUE)
scca_cca_tr_cor_se <- rowSds(simplify2array(scca_cca_tr_ro['cors', ]), na.rm = TRUE) / sqrt(dim(scca_cca_tr_ro)[2])

scca_corr <- data.frame(
  mode = as.factor(1:nMode),
  rval = scca_cca_tr_cor,
  se   = scca_cca_tr_cor_se
)
scca_corr_reorder <- scca_corr %>% arrange(desc(rval)) %>% mutate(neworder = 1:nMode)
scca_corr_reorder


# Permutation Analysis ---------------------------------------------------------
# run sCCA with permuted externalizing symptoms
nPerms <- 5000
set.seed(1111)
ext_resid_perm_tr <- rlply(nPerms, data_train$ext[sample(nrow(data_train$ext)), ])
scca.perm.cca <- sapply(
  ext_resid_perm_tr, function(y_perm) {
    out <- ccaDWpermorder(data_train$fconn, y_perm, pen_x_selected, pen_z_selected, nMode, scca_cand_tr)
  }
)
scca_perm_tr_rval <- simplify2array(scca.perm.cca["cors", ])
scca_perm_tr_pval <- sapply(
  seq_along(scca_corr$rval), function(x) {
    sum(scca_perm_tr_rval[x, ] >= scca_corr$rval[x], na.rm = TRUE) / sum(!is.na(scca_perm_tr_rval[x, ]))
  }
)
scca_perm_tr_pval_adj <- p.adjust(scca_perm_tr_pval, method = "bonferroni")
scca_perm_tr_pass <- which(scca_perm_tr_pval_adj < 0.05)

dim.match <- sapply(
  seq_along(1:length(scca_perm_tr_pass)), function(x) {
    which(scca_perm_tr_pass == scca_corr_reorder$mode[x])
  }
)
# add second component which survives uncorr


# plot and save correlation plot
peak <- ceiling(scca_corr_reorder$rval[1] * 10) / 10
label1 <- sprintf("italic(p[bonf]) == %0.3f", scca_perm_tr_pval_adj[scca_perm_tr_pass[dim.match[1]]])
label2 <- sprintf("italic(p[bonf]) == %0.3f", scca_perm_tr_pval_adj[scca_perm_tr_pass[dim.match[2]]])
corrplot <- scca_corr_reorder %>% 
  mutate(ymax = .$rval + .$se) %>%
  mutate(ymin = .$rval - .$se) %>%
  ggplot(aes(neworder, rval)) +
  geom_bar(width = 0.75, stat = "identity", fill = c("gray20", "gray20", "gray60")) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.3) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, peak), breaks=seq(0, peak, 0.05)) +
  labs(x = "", y = "CCA Correlation") +
  theme_classic(base_size = 12) +
  annotate("text", x = 1, y = scca_corr_reorder$rval[1] + 0.02, label = label1, parse = TRUE, size = 4) +
  annotate("text", x = 2, y = scca_corr_reorder$rval[2] + 0.02, label = label2, parse = TRUE, size = 4) +
  coord_cartesian(ylim=c(0.2, peak))
corrplot
ggsave(file.path("figures", "adhd_train_corrplot.png"), corrplot, width = 4, height = 3, dpi = 400)


# Bootstrap Analysis for Stability ---------------------------------------------
scca_boot = list(NULL)
nBoots = 5000
for (ith_boot in 1:nBoots) {
  set.seed(ith_boot)
  bootid <- createResample(data_tr$Age, list = FALSE, times = 1)
  fconn_boot <- data_train$fconn[bootid, ]
  ext_boot   <- data_train$ext[bootid, ]
  scca_boot[[ith_boot]]  <- run_scca(fconn_boot, ext_boot, pen_x_selected, pen_z_selected, nMode)
}
scca_boot.ro <- lapply(1:nBoots, function(i) reorderCCA(scca_boot[[i]], scca_cand_tr, nMode))
scca_boot.u  <- lapply(scca_perm_tr_pass, function(x) sapply(1:nBoots, function(i) scca_boot.ro[[i]]$u[, x])) #extract loadings on connectivity features
scca_boot.v  <- lapply(scca_perm_tr_pass, function(x) sapply(1:nBoots, function(i) scca_boot.ro[[i]]$v[, x])) #extract loadings on clinical features
u.boot.plot  <- lapply(seq_along(scca_perm_tr_pass), function(x) bootplot_u(scca_cand_tr$u[, scca_perm_tr_pass[x]], scca_boot.u[[x]] ))  #apply 99.5% confidence interval
#v.boot.plot  <- lapply(seq_along(scca_perm_tr_pass), function(x) bootplot(scca_cand_tr$v[,scca_perm_tr_pass[x]], scca_boot.v[[x]] )) #apply 95% confidence interval


# Visualization: Beh -----------------------------------------------------------
item_labels <- read_csv(file.path(path_wd, "data", "raw", "pnc_ext_label.csv")) %>% pull("items")
a <- bootplot4beh(scca_cand_tr$v[, scca_perm_tr_pass[dim.match[1]]], scca_boot.v[[1]], item_labels)
b <- bootplot4beh(scca_cand_tr$v[, scca_perm_tr_pass[dim.match[2]]], scca_boot.v[[2]], item_labels)
#v_plot <- ggarrange(a$plot + rremove("x.text"), b$plot + rremove("x.text"), ncol = 1, nrow = 2) %>%
#  ggexport(
#    filename = file.path("figures", "adhd_beh_features.pdf"),
#    width    = 5,
#    height   = 4
#  )

ggsave(file.path("figures", "adhd_mode1_behavior.png"), plot = a$plot + rremove("y.text"), height = 2, width = 3)
ggsave(file.path("figures", "adhd_mode2_behavior.png"), plot = b$plot + rremove("y.text"), height = 2, width = 3)


# Visulization: FC -------------------------------------------------------------
# Mode 1 FC Visualization ------------------------------------------------------
fconn_feature <- matrix(0, 268, 3)
fconn_feature[, 1] <- u.boot.plot[[dim.match[1]]]$load_nm[1:268 + 268 * 0] %>% abs()
fconn_feature[, 2] <- u.boot.plot[[dim.match[1]]]$load_nm[1:268 + 268 * 1] %>% abs()
fconn_feature[, 3] <- u.boot.plot[[dim.match[1]]]$load_nm[1:268 + 268 * 2] %>% abs()
shen_labels <- read_excel("Shen_idx.xls")

labels <- 1:268
components <- c("MCC", "pre-SMA", "Striatum")
rownames(fconn_feature) <- labels
colnames(fconn_feature) <- components

#tmp <- u.boot.plot[[dim.match[1]]]$load_nm %>% abs() %>% sort(decreasing = TRUE)  
#cutoff <- tmp[round(dim(data_tr)[2] * 0.05)]
#fconn_feature[fconn_feature < cutoff] <- 0

network <- table(shen_labels$YeoWithSubc)
network_name_org <- rownames(network)
network_name <- c("VIS", "SMN", "DAN", "VAN", "LIM", "FPN", "DMN", "SBN", "CBN")

sumNet <- matrix(0, 9, 3)
for (ith_net in 1:9) {
  for (ith_roi in 1:3){
    index <- which(shen_labels$YeoWithSubc == ith_net)
    sumNet[ith_net, ith_roi] <- fconn_feature[index, ith_roi] %>% sum() / network[ith_net]
  }
}

factors <- c(components, "Empty1", network_name, "Empty2")
n_factors = length(factors)
# color setting
pal1 <- brewer.pal(n = 10, name = 'Paired')
pal2 <- brewer.pal(n = 8, name = 'Dark2')
combined_pal <- c(pal1[c(2, 3, 1, 4)], pal2[c(2, 1, 3, 6, 8)])
color <- c("black", "black", "black", "white", combined_pal, "white")
colpal <- c("darkseagreen4", "goldenrod3", "deepskyblue4")


# Mode 2 FC Visualization ------------------------------------------------------
fconn_feature <- matrix(0, 268, 3)
fconn_feature[, 1] <- u.boot.plot[[dim.match[2]]]$load_nm[1:268 + 268 * 0] %>% abs()
fconn_feature[, 2] <- u.boot.plot[[dim.match[2]]]$load_nm[1:268 + 268 * 1] %>% abs()
fconn_feature[, 3] <- u.boot.plot[[dim.match[2]]]$load_nm[1:268 + 268 * 2] %>% abs()
shen_labels <- read_excel("Shen_idx.xls")

labels <- 1:268
components <- c("MCC", "pre-SMA", "Striatum")
rownames(fconn_feature) <- labels
colnames(fconn_feature) <- components

#tmp <- u.boot.plot[[dim.match[1]]]$load_nm %>% abs() %>% sort(decreasing = TRUE)  
#cutoff <- tmp[round(dim(data_tr)[2] * 0.05)]
#fconn_feature[fconn_feature < cutoff] <- 0

sumNet <- matrix(0, 9, 3)
for (ith_net in 1:9) {
  for (ith_roi in 1:3){
    index <- which(shen_labels$YeoWithSubc == ith_net)
    sumNet[ith_net, ith_roi] <- fconn_feature[index, ith_roi] %>% sum() / network[ith_net]
  }
  
}

coul <- brewer.pal(3, "Set1")
colors_border <- coul
colors_in <- alpha(colpal, 0.3)
sumNet2 <- cbind(rep(0.4, 9), rep(0, 9), sumNet)
rownames(sumNet2) <- network_name

pdf(file.path("figures", "adhd_network_mode2_features.pdf"), width = 6, height = 6)
sumNet2 %>%
  t() %>% 
  as.data.frame() %>%
  radarchart(
    axistype = 1,
    pcol = colpal, pfcol= alpha(colpal, 0.3), plwd = 2.5, plty = 1,
    #custom the grid
    cglcol = "grey", cglty = 1.2, axislabcol = "grey", caxislabels = seq(0, 0.4, 0.1), cglwd = 0.5,
    #custom labels
    vlcex = 1
  )
dev.off()



network <- table(shen_labels$Network_canonical)
network_name_org <- rownames(network)

lobes <- table(shen_labels$`Lobe No`)
lobe_labels <- c("Prefrontal", "MotorStrip", "Insula", "Parietal", "Temporal",
                 "Occipital", "Limbic", "Cerebellum", "Subcortical", "Brainstem")

factors <- c(components, "Empty1", network_name, "Empty2")
n_factors = length(factors)
# color setting
pal1 <- brewer.pal(n = 10, name = 'Paired')
pal2 <- brewer.pal(n = 8, name = 'Dark2')
combined_pal <- c(pal1[c(2, 3, 1)], pal2[c(2, 1, 3, 6, 8)])
lobe_pal <- c(pal1[c(2, 3, 1, 4, 5)], pal2[c(2, 1, 3, 6, 8)])
color <- c("black", "black", "black", "white", combined_pal, "white")
colpal <- c("darkseagreen4", "goldenrod3", "deepskyblue4")







# sCCA on Replication ----------------------------------------------------------
scca_main_te <- run_scca(data_test$fconn, data_test$ext, pen_x_selected, pen_z_selected)

EVs_te <- scca_main_te$d ^ 2 / sum(scca_main_te$d ^ 2)
EVs_reorder_te <- sort(EVs_te, decreasing = TRUE, index.return = TRUE)
EVs_df_te <- data.frame(Mode = as.factor(1:nItems), EV = EVs_reorder_te$x)
maxEV_te <- ceiling(max(EVs_te) * 100) / 100 + 0.005

EVs_screeplot_te <- EVs_df_te %>%
  mutate(highlight_flag = ifelse(Mode %in% 1:3, TRUE, FALSE)) %>% 
  ggplot(aes(Mode, EV, color = highlight_flag)) +
  geom_point(stat = 'identity', aes(size = EV)) +
  scale_color_manual(values = c("grey80", "grey35")) +
  geom_hline(yintercept = expectedEV, linetype = "dashed", color = "grey50") +
  scale_x_discrete(breaks = seq(1, 22, 1)) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, maxEV_te),
    breaks = seq(0, maxEV_te, 0.02)
  ) +
  labs(x = "Mode", y = "Explained Variance") +
  theme_classic(base_size = 12) +
  theme(legend.position = 'none')
EVs_screeplot_te
ggsave(file.path("figures", "adhd_replication_screeplot.png"), EVs_screeplot_te, width = 6, height = 4, dpi = 400)

# run sCCA with permutation
scca_cand_te <- run_scca(data_test$fconn, data_test$ext, pen_x_selected, pen_z_selected, nMode)

# create 3 folds 100 times in the sub-training set to obtain average
set.seed(1111)
resample_idx_te     <- createDataPartition(data_te$EXT, p = 0.667, list = TRUE, times = 100)
fconn_10resample_te <- mclapply(resample_idx_te, function(x) data_test$fconn[x,])
ext_10resample_te   <- mclapply(resample_idx_te, function(x) data_test$ext[x,])

scca_cca_te  <- mclapply(
  seq_along(resample_idx_te), function(i) {
    run_scca(fconn_10resample_te[[i]], ext_10resample_te[[i]], pen_x_selected, pen_z_selected, nMode)
  }
)
scca_cca_te_ro     <- sapply(scca_cca_te, function(x) reorderCCA(x, scca_cand_te, nMode))
scca_cca_te_cor    <- rowMeans(simplify2array(scca_cca_te_ro['cors', ]), na.rm = TRUE)
scca_cca_te_cor_se <- rowSds(simplify2array(scca_cca_te_ro['cors', ]), na.rm = TRUE) / sqrt(dim(scca_cca_te_ro)[2])

scca_corr_te <- data.frame(
  mode = as.factor(1:nMode),
  rval = scca_cca_te_cor,
  se   = scca_cca_te_cor_se
)
scca_corr_te_reorder <- scca_corr_te %>% arrange(desc(rval)) %>% mutate(neworder = 1:nMode)
scca_corr_te_reorder


# Permutation Analysis ---------------------------------------------------------
# run sCCA with permuted externalizing symptoms
set.seed(1111)
ext_resid_perm_te <- rlply(nPerms, data_test$ext[sample(nrow(data_test$ext)), ])
scca.perm.cca_test <- sapply(
  ext_resid_perm_te, function(y_perm) {
    out <- ccaDWpermorder(data_test$fconn, y_perm, pen_x_selected, pen_z_selected, nMode, scca_cand_te)
  }
)
scca_perm_te_rval <- simplify2array(scca.perm.cca_test["cors", ])
scca_perm_te_pval <- sapply(
  seq_along(scca_corr_te$rval), function(x) {
    sum(scca_perm_te_rval[x, ] >= scca_corr_te$rval[x], na.rm = TRUE) / sum(!is.na(scca_perm_te_rval[x, ]))
  }
)
scca_perm_te_pval_adj <- p.adjust(scca_perm_te_pval, method = "bonferroni")
scca_perm_te_pass <- which(scca_perm_te_pval_adj < 0.05)
scca_perm_te_pass <- c(2, 3)

dim.match_te <- sapply(
  seq_along(1:length(scca_perm_te_pass)), function(x) {
    which(scca_perm_te_pass == scca_corr_te_reorder$mode[x])
  }
)


# plot text correlation bar plot
peak <- ceiling(scca_corr_te_reorder$rval[1] * 10) / 10 + 0.05
label1 <- sprintf("italic(p[bonf]) == %0.3f", scca_perm_te_pval_adj[scca_perm_te_pass[dim.match_te[1]]])
label2 <- sprintf("italic(p[bonf]) == %0.3f", scca_perm_te_pval_adj[2])
corrplot_te <- scca_corr_te_reorder %>% 
  mutate(ymax = .$rval + .$se) %>%
  mutate(ymin = .$rval - .$se) %>%
  ggplot(aes(neworder, rval)) +
  geom_bar(width = 0.75, stat = "identity", fill = c("gray20", "gray60", "gray60")) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.3) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, peak), breaks=seq(0, peak, 0.05)) +
  labs(x = "", y = "CCA Correlation") +
  theme_classic(base_size = 12) +
  annotate("text", x = 1, y = scca_corr_te_reorder$rval[1] + 0.02, label = label1, parse = TRUE, size = 4) +
  annotate("text", x = 2, y = scca_corr_te_reorder$rval[2] + 0.02, label = label2, parse = TRUE, size = 4) +
  #annotate("text", x = 3, y = scca_corr_te_reorder$rval[3] + 0.02, label = label3, parse = TRUE, size = 4) +
  coord_cartesian(ylim=c(0.35, peak))
corrplot_te
ggsave(file.path("figures", "adhd_replication_corrplot.pdf"), corrplot_te, width = 4, height = 3)

# Bootstrap Analysis for Stability ---------------------------------------------
scca_boot_te = list(NULL)
for (ith_boot in 1:nBoots) {
  set.seed(ith_boot)
  bootid_te <- createResample(data_te$Age, list = FALSE, times = 1)
  fconn_boot_te <- data_test$fconn[bootid_te, ]
  ext_boot_te   <- data_test$ext[bootid_te, ]
  scca_boot_te[[ith_boot]]  <- run_scca(fconn_boot_te, ext_boot_te, pen_x_selected, pen_z_selected, nMode)
}
scca_boot_te.ro <- lapply(1:nBoots, function(i) reorderCCA(scca_boot_te[[i]], scca_cand_te, nMode))
scca_boot_te.u  <- lapply(scca_perm_te_pass, function(x) sapply(1:nBoots, function(i) scca_boot_te.ro[[i]]$u[, x])) #extract loadings on connectivity features
scca_boot_te.v  <- lapply(scca_perm_te_pass, function(x) sapply(1:nBoots, function(i) scca_boot_te.ro[[i]]$v[, x])) #extract loadings on clinical features
u.boot_te.plot  <- lapply(seq_along(scca_perm_te_pass), function(x) bootplot_u(scca_cand_te$u[, scca_perm_te_pass[x]], scca_boot_te.u[[x]] ))  #apply 99.5% confidence interval
v.boot_te.plot  <- lapply(seq_along(scca_perm_te_pass), function(x) bootplot(scca_cand_te$v[,scca_perm_te_pass[x]], scca_boot_te.v[[x]] )) #apply 95% confidence interval

date <- gsub("-", "", Sys.Date())
save.image(file = sprintf("%s_adhd_analysis.RData", date))





