HexBin_Signal <- function(data, pixel_info) {

  hexbin_list <- list()
  
  for (scan in unique(pixel_info$run)) {
    
    for (hexbin in unique(pixel_info$bin)) {
      idx_run <- pixel_info$run == scan
      idx_hexbin <- pixel_info$bin == hexbin
      idx_on_tissue <- pixel_info$on_tissue == T
      
      hexbin_idx <- Reduce("&", list(idx_run, idx_hexbin, idx_on_tissue))
      
      hexbin_df <- data[, hexbin_idx, drop = F]
      
      if (ncol(hexbin_df) == 0) next
      mean_vec <- rowMeans(hexbin_df)
      
      temp <- data.frame(mean_intensity = mean_vec, row.names = names(mean_vec))
      colnames(temp) <- paste(scan, hexbin, sep = "@")
      hexbin_list[[paste(scan, hexbin, sep = "@")]] <- temp
    }
    
    
  }
  
  hexbin_df_all <- do.call(cbind, hexbin_list)
  
  return(hexbin_df_all)
}