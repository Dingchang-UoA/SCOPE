calculate_weights <- function(j, coords, span) {
  
  # Step 1: calculate distances
  dists <- sqrt((coords$x[j] - coords$x)^2 + (coords$y[j] - coords$y)^2)
  
  # Step 2: raw weights using triangle kernel
  raw_weights <- pmax(0, 1 - dists / span)
  
  # Step 3: normalize the weights (important for correct statistics)
  normalized_weights <- raw_weights / sum(raw_weights)
  
  return(normalized_weights)[j]
}
