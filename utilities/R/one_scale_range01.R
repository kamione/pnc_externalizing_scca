one_scale_range01 <- function(input) {
  max_val <- max(input)
  min_val <- min(input)
  output <- (input - min_val) / (max_val - min_val)
  return(output)
}