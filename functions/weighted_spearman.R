weighted_spearman <- function(x, y, weights) {
  rx <- rank(x)
  ry <- rank(y)
  
  w_mean_rx <- sum(weights * rx)
  w_mean_ry <- sum(weights * ry)
  
  num <- sum(weights * rx * ry) - w_mean_rx * w_mean_ry
  
  denom_x <- sum(weights * rx^2) - w_mean_rx^2
  denom_y <- sum(weights * ry^2) - w_mean_ry^2
  
  denom <- sqrt(denom_x * denom_y)
  
  if (denom == 0) return(NA)  # to avoid division by zero
  
  corr <- num / denom
  return(corr)
}