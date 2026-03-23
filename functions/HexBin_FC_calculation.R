HexBin_FC_calculation <- function(data = hexbin_average_flt,
                                  group1,
                                  group2
                                    ){
  hbn_mean_raw <- hexbin_average_flt %>% select(matches(paste0("^", group1)))
  hbn_mean_raw <- hbn_mean_raw+1
  hbo_mean_raw <- hexbin_average_flt %>% select(matches(paste0("^", group2)))
  hbo_mean_raw <- hbo_mean_raw+1
  
  hbn_mean_info <- as.data.frame(str_split_fixed(colnames(hbn_mean_raw), "@", 2))
  colnames(hbn_mean_info) <- c("run", "bin")
  hbo_mean_info <- as.data.frame(str_split_fixed(colnames(hbo_mean_raw), "@", 2))
  colnames(hbo_mean_info) <- c("run", "bin")
  
  fc_bin_raw <- data.frame(mz = rownames(hbn_mean_raw))
  
  for (bin_idx in unique(hbo_mean_info$bin)) {
    idx1 <- which(hbo_mean_info$bin == bin_idx)
    idx2 <- which(hbn_mean_info$bin == bin_idx)
    
    
    # If both idx1 and idx2 exist, calculate fc_bin
    if (length(idx1) > 0 && length(idx2) > 0) {
      
      fc_bin_raw[, paste0(bin_idx)] <- hbo_mean_raw[, idx1] / hbn_mean_raw[, idx2]
    }
  }
  rownames(fc_bin_raw) <- fc_bin_raw$mz
  fc_bin_raw_df <- fc_bin_raw[,-1]
  return(fc_bin_raw_df)
}