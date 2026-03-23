Average_HexBin_intensity <- function(data) {
  info <- as.data.frame(str_split_fixed(colnames(data), "@", 2))
  colnames(info) <- c("run", "bin")
  info$treatment <-  sub("_.*", "", info$run)
  
  hexbin_average <- list()
  
  for (m in unique(info$treatment)) {
    for (n in unique(info$bin)) {
      idx_treatment <- info$treatment == m
      idx_hexbin <- info$bin == n 
      
      idx <- '&'(idx_treatment, idx_hexbin)
      hexbin_df <- data[, idx, drop = F]
      
      if (ncol(hexbin_df) == 0) next
      mean_vec <- rowMeans(hexbin_df)
      
      temp <- data.frame(mean_intensity = mean_vec, row.names = names(mean_vec))
      colnames(temp) <- paste(m, n, sep = "@")
      hexbin_average[[paste(m, n, sep = "@")]] <- temp
    }
  }
  
  
  hexbin_average_df <- do.call(cbind, hexbin_average)
  return(hexbin_average_df)
}

