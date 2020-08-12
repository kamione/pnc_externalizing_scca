one_screeplot <- function(parallel, main.title = NULL) {
  # Info
  # Function:
  #   
  # Arguements:
  #   parallel: output of fa.parallel function
  #
  # Output
  #   ggplot figure
  
  # creates a data frame from observed eigenvalue data
  if (is.null(main.title) == 1) {
    main.title = "Parallel Analysis Scree Plot" 
  }
  
  obs.dat <- data.frame(parallel$fa.values)
  obs.dat$type <- c("Observed Data")
  obs.dat$num  <- as.numeric(c(row.names(obs.dat)))
  colnames(obs.dat) <- c("eigenvalue", "type", "num")
  
  # obtains the 95th percentile of the simulated eigenvalues
  percentile <- apply(parallel$values, 2, function(x) quantile(x,.95))
  min <- as.numeric(nrow(obs.dat))
  min <- (4 * min) - (min - 1)
  max <- as.numeric(nrow(obs.dat))
  max <- 4 * max
  percentile1 <- percentile[min:max]
  
  # creates a data frame from simulated eigenvalue data
  sim.dat <- data.frame(percentile1)
  sim.dat$type <- c("Simulated Data (95th percentile)")
  sim.dat$num  <- as.numeric(c(row.names(obs.dat)))
  colnames(sim.dat) <- c("eigenvalue", "type", "num")
  
  # combines the observed and simulated eigenvalue data
  eigen.dat <- rbind(obs.dat, sim.dat)
  
  # creates a theme according to APA format
  apatheme <- theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 16, face = "bold",
                                    hjust = 0.5, vjust = 1),
          text = element_text(family = "Arial"),
          legend.title = element_blank(),
          legend.position = c(.75, .88),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"))
  
  # plots the scree plot with observed and simulated data
  p <- ggplot(eigen.dat, aes(x = num, y = eigenvalue, shape = type)) +
    geom_hline(yintercept = 1) +
    geom_vline(xintercept = parallel$nfact, linetype = "dashed") +
    geom_line() +
    geom_point(size = 3) +
    ggtitle(main.title) +
    scale_y_continuous(name = "Eigenvalue",
                       breaks = seq(0, max(eigen.dat$eigenvalue), 2)) +
    scale_x_continuous(name = "Factor Number",
                       breaks = min(eigen.dat$num):max(eigen.dat$num)) +
    scale_shape_manual(values = c(16, 1)) +
    apatheme
  return(p)
}