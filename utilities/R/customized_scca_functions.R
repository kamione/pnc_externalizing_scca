# Information ------------------------------------------------------------------
# Version    : 0.1.0
# Date       : 13-April-2020
# Description: 

# Load Packages ----------------------------------------------------------------
library(tidyverse)
library(PMA)

# Functions --------------------------------------------------------------------
run_scca <- function(X, Z, pen_x, pen_z, nMode = NULL) {
  nX <- dim(X)[2]
  nZ <- dim(Z)[2]
  if ( is.null(nMode) ) {
    nMode <- min(c(nX, nZ))
  }
  output <- PMA::CCA(
    x        = X,
    z        = Z,
    typex    = "standard",
    typez    = "standard",
    penaltyx = pen_x,
    penaltyz = pen_z,
    K        = nMode,
    niter    = 100,
    trace    = FALSE
  )
  return(output)
}

searh_hyperparameters <- function(X, Z, nPerms = 100) {
  output <- PMA::CCA.permute(
    x         = X,
    z         = Z,
    typex     = "standard",
    typez     = "standard",
    penaltyxs = seq(0.1, 1, 0.05),
    penaltyzs = seq(0.1, 1, 0.05),
    nperms    = nPerms,
    niter     = 20,
    trace     = FALSE
  )
  return(output)
}

reorderCCA_new <- function(res, org, k) {
  u <- org$u
  v <- org$v
  
  res.match<-apply(abs(cor(v,res$v[,1:k])),1,which.max)
  res.cor <- apply(abs(cor(v,res$v[,1:k])),1,max)
  res.match.count <- length(unique(res.match))
  
  u.od <- res$u[,res.match]
  v.od <- res$v[,res.match]
  cors.final <- res$cors[res.match]
  
  org.sign <- sign(colMeans(sign(v)))
  res.sign <- sign(colMeans(sign(v.od)))
  
  sign.prod <- org.sign *res.sign
  
  u.final <- t(t(u.od) * sign.prod)
  v.final <- t(t(v.od) * sign.prod)
  
  res.one.reorder <- list(u= u.final ,v= v.final, cors = cors.final, pos = res.match, dimcor = res.cor)
  
  out <- res.one.reorder
}

NullSpace <- function (A) {
  m <- dim(A)[1]; n <- dim(A)[2]
  ## QR factorization and rank detection
  QR <- base::qr.default(A)
  r <- QR$rank
  ## cases 2 to 4
  if ((r < min(m, n)) || (m < n)) {
    R <- QR$qr[1:r, , drop = FALSE]
    P <- QR$pivot
    F <- R[, (r + 1):n, drop = FALSE]
    I <- base::diag(1, n - r)
    B <- -1.0 * base::backsolve(R, F, r)
    Y <- base::rbind(B, I)
    X <- Y[base::order(P), , drop = FALSE]
    return(X)
  }
  ## case 1
  return(base::matrix(0, n, 1))
}

bootstats2 <- function(org,bootdata, uplim,lwlim){
  boot.stats <- sapply(seq_along(1:length(org)), function(i) org[i] - quantile(bootdata[i,] - org[i],c(uplim,lwlim),na.rm =T))
  boot.stats <- as.data.frame(t(boot.stats))
  colnames(boot.stats)<-c("low","high")
  boot.stats$ci <- boot.stats$high - boot.stats$low
  boot.stats$load <- org
  boot.stats$fea <- 1:length(org)
  boot.stats
}

bootplot4beh <- function(org, boot, labels) {
  btst <- bootstats2(org, boot, 0.975, 0.025)
  
  #p <- btst %>%
  #  ggplot(aes(fea, load)) +
  #    geom_point(aes(color = low * high > 0), size = 4) +
  #    scale_color_manual(values = c("darkgrey", "tomato2")) +
  #    geom_errorbar(aes(ymax = high, ymin = low),  width = 0.25) +
  #    labs(x = "", y = "Loading") +
  #    scale_x_discrete(limits = 1:22, labels = labels) +
  #    theme_classic() +
  #    theme(legend.position = "none") +
  #    theme(axis.text.x = element_text(
  #      hjust = 0.95,
  #      vjust = 1,
  #      face  = "bold",
  #      angle = 60,
  #      color = c(rep("#3886aa", 9), rep("#ff3609", 8), rep("#ffc109", 5)))
  #    )
  type = c(rep(1, 9), rep(2, 8), rep(3, 5))
  color = c("#3886aa", "#ff3609", "#ffc109")
  tmp <- btst %>%
    mutate(diagnosis = type) %>% 
    filter(.$low * .$high > 0) %>%
    mutate(label = factor(labels[fea], levels = rev(labels[fea])))
  
  p <- tmp %>%
    #ggplot(aes(x = label, y = abs(load), color = as.factor(diagnosis))) +
    ggplot(aes(x = label, y = abs(load)))+
      geom_point(size = 5) +
      #scale_color_manual(values=color[unique(tmp$diagnosis)]) +
      geom_errorbar(aes(ymax = abs(high), ymin = abs(low)),  width = 0.2) +
      scale_y_continuous(limits = c(0, 1.5), breaks=seq(0, 1.5, by = 0.2)) +
      labs(x = "", y = "Canonical Loading") +
      coord_flip() +
      theme_classic() +
      theme(legend.position = "none")
  
  fea <- which(btst$low * btst$high > 0)
  load <- org[fea]
  fea <- fea[order(-abs(load))]
  load <- load[order(-abs(load))]
  load_nm <- org / btst$ci
  load_nm_fea <- load_nm[fea]
  out <- list(plot = p, fea = fea, load = load, ci = btst$ci[fea], high = btst$high[fea], low = btst$low[fea],load_nm = load_nm, load_nm_fea = load_nm_fea)
  return(out)
}

plot_fconn_3seeds <- function() {
  # Visulization: FC -------------------------------------------------------------
  fconn_feature <- matrix(0, 268, nSeeds)
  fconn_feature[, 1] <- u.boot.plot[[2]]$load_nm[1:268]
  fconn_feature[, 2] <- u.boot.plot[[2]]$load_nm[1:268 + 268]
  shen_labels <- read_excel("Shen_idx.xls")
  
  labels <- 1:268
  components <- c("Amygdala", "IPL")
  rownames(fconn_feature) <- labels
  colnames(fconn_feature) <- components
  
  network <- table(shen_labels$Network_canonical)
  network_name_org <- rownames(network)
  network_name <- c("DMN", "FPN", "MF", "MOT", "SBC", "VAS", "VI", "VII")
  
  factors <- c(components, "Empty1", network_name, "Empty2")
  n_factors = length(factors)
  # color setting
  pal1 <- brewer.pal(n = 10, name = 'Paired')
  pal2 <- brewer.pal(n = 8, name = 'Dark2')
  combined_pal <- c(pal1[c(2, 3, 1)], pal2[c(2, 1, 3, 6, 8)])
  color <- c("black", "black", "white", combined_pal, "white")
  
  edgecollist <- c("gray15", "goldenrod", "navy")
  
  
  # Circular Plot ----------------------------------------------------------------
  pdf(file.path("figures", "dbd_05_fconn_features.pdf"), width = 4, height = 4)
  
  par(mar = c(0, 0, 0, 0))
  circos.par(start.degree = 275)
  mat_layout = cbind(rep(0, n_factors), c(4, 4, 3, as.vector(network) / 8, 3))
  circos.initialize(factors, xlim = mat_layout)
  # create empty outer ring for the highlight
  circos.trackPlotRegion(
    track.index = 1,
    ylim = c(0, 1), 
    track.height = 0.3, 
    bg.col = "white",
    panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, 
                  facing = "bending", niceFacing = TRUE, cex = 0.4, col = "white")
    },
    bg.border = "white"
  )
  
  # second ring
  circos.trackPlotRegion(
    track.index = 2,
    ylim = c(0, 1), 
    track.height = 0.07, 
    bg.col = color, 
    panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, 
                  facing = "bending", niceFacing = TRUE, cex = 0.7, col = "white",
                  font = 2)
    },
    bg.border = "white"
  )
  # links from factor a to b
  node_flag <- matrix(0, 268, 1)
  max_pval <- max(fconn_feature)
  for (sign in 1:2) {
    edgecol = edgecollist[sign]
    idx_flag <- matrix(0, 8, 1) # create a matrix to store flag of networks
    for (ith_node in 1:268) {
      netname <- shen_labels$Network_canonical[ith_node]
      idx <- match(netname, network_name_org)
      idx_flag[idx] <- idx_flag[idx] + 1
      if (fconn_feature[ith_node, sign] == 0) { 
        next
      } else {
        node_flag[ith_node] <- node_flag[ith_node] + 1
        alpha_value <- as.numeric(fconn_feature[ith_node, sign]) / max_pval * 0.9
      }
      edge_pal <- adjustcolor(edgecol, alpha.f = alpha_value)
      startend <- c((idx_flag[idx] - 1), idx_flag[idx]) / 8
      circos.link(factors[sign], c(0, 4), factors[(idx + 3)], startend, col = edge_pal, border = NA)
    }
  }
  circos.clear()
  dev.off()
  
  # mode 2
  fconn_feature <- matrix(0, 268, 2)
  fconn_feature[, 1] <- u.boot.plot[[1]]$load_nm[1:268]
  fconn_feature[, 2] <- u.boot.plot[[1]]$load_nm[1:268 + 268]
  shen_labels <- read_excel("Shen_idx.xls")
  
  labels <- 1:268
  components <- c("Amygdala", "IPL")
  rownames(fconn_feature) <- labels
  colnames(fconn_feature) <- components
  
  network <- table(shen_labels$Network_canonical)
  network_name_org <- rownames(network)
  network_name <- c("DMN", "FPN", "MF", "MOT", "SBC", "VAS", "VI", "VII")
  
  factors <- c(components, "Empty1", network_name, "Empty2")
  n_factors = length(factors)
  # color setting
  pal1 <- brewer.pal(n = 10, name = 'Paired')
  pal2 <- brewer.pal(n = 8, name = 'Dark2')
  combined_pal <- c(pal1[c(2, 3, 1)], pal2[c(2, 1, 3, 6, 8)])
  color <- c("black", "black", "white", combined_pal, "white")
  
  edgecollist <- c("gray15", "goldenrod", "navy")
  
  pdf(file.path("figures", "dbd_05_mode2_fconn_features.pdf"), width = 4, height = 4)
  
  par(mar = c(0, 0, 0, 0))
  circos.par(start.degree = 275)
  mat_layout = cbind(rep(0, n_factors), c(4, 4, 3, as.vector(network) / 8, 3))
  circos.initialize(factors, xlim = mat_layout)
  # create empty outer ring for the highlight
  circos.trackPlotRegion(
    track.index = 1,
    ylim = c(0, 1), 
    track.height = 0.3, 
    bg.col = "white",
    panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, 
                  facing = "bending", niceFacing = TRUE, cex = 0.4, col = "white")
    },
    bg.border = "white"
  )
  
  # second ring
  circos.trackPlotRegion(
    track.index = 2,
    ylim = c(0, 1), 
    track.height = 0.07, 
    bg.col = color, 
    panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, 
                  facing = "bending", niceFacing = TRUE, cex = 0.7, col = "white",
                  font = 2)
    },
    bg.border = "white"
  )
  # links from factor a to b
  node_flag <- matrix(0, 268, 1)
  max_pval <- max(fconn_feature)
  for (sign in 1:2) {
    edgecol = edgecollist[sign]
    idx_flag <- matrix(0, 8, 1) # create a matrix to store flag of networks
    for (ith_node in 1:268) {
      netname <- shen_labels$Network_canonical[ith_node]
      idx <- match(netname, network_name_org)
      idx_flag[idx] <- idx_flag[idx] + 1
      if (fconn_feature[ith_node, sign] == 0) { 
        next
      } else {
        node_flag[ith_node] <- node_flag[ith_node] + 1
        alpha_value <- as.numeric(fconn_feature[ith_node, sign]) / max_pval * 0.9
      }
      edge_pal <- adjustcolor(edgecol, alpha.f = alpha_value)
      startend <- c((idx_flag[idx] - 1), idx_flag[idx]) / 8
      circos.link(factors[sign], c(0, 4), factors[(idx + 3)], startend, col = edge_pal, border = NA)
    }
  }
  circos.clear()
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}




