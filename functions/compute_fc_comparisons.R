compute_fc_comparisons <- function(
    hexbin_average_df,
    comparisons,
    standard_lens,
    target_mz = 306.0773
) {
  

  # Extract column info
  info <- as.data.frame(stringr::str_split_fixed(colnames(hexbin_average_df), "@", 2))
  colnames(info) <- c("treatment", "bin")
  
  comparisons_lt <- vector("list", length(comparisons))
  
  for (i in seq_along(comparisons)) {
    
    comp <- comparisons[[i]]
    
    # Filter relevant columns
    hexbin_average_flt <- hexbin_average_df[, info$treatment %in% comp, drop = FALSE]
    
    # Split into two groups
    grp1_raw <- hexbin_average_flt %>% dplyr::select(dplyr::matches(paste0("^", comp[1])))
    grp2_raw <- hexbin_average_flt %>% dplyr::select(dplyr::matches(paste0("^", comp[2])))
    
    # Extract bin info
    grp1_info <- as.data.frame(stringr::str_split_fixed(colnames(grp1_raw), "@", 2))
    colnames(grp1_info) <- c("run", "bin")
    
    grp2_info <- as.data.frame(stringr::str_split_fixed(colnames(grp2_raw), "@", 2))
    colnames(grp2_info) <- c("run", "bin")
    
    # Initialize FC table
    fc_bin <- data.frame(mz = rownames(grp1_raw))
    
    # Compute fold change per bin
    for (bin_idx in unique(grp2_info$bin)) {
      idx1 <- which(grp2_info$bin == bin_idx)
      idx2 <- which(grp1_info$bin == bin_idx)
      
      if (length(idx1) > 0 && length(idx2) > 0) {
        fc_bin[[as.character(bin_idx)]] <- grp2_raw[, idx1] / grp1_raw[, idx2]
      }
    }
    
    # Clean + format
    rownames(fc_bin) <- fc_bin$mz
    fc_mat <- as.matrix(fc_bin[, -1, drop = FALSE])
    
    # Extract target m/z
    fc_target <- as.data.frame(t(fc_mat[which(rownames(fc_mat) == target_mz), , drop = FALSE]))
    colnames(fc_target) <- "FC"
    fc_target$bin <- as.numeric(rownames(fc_target))
    fc_target$comparison <- paste0(comp[1], " vs ", comp[2])
    
    # Join to spatial map
    standard_lens_hex <- standard_lens %>%
      dplyr::left_join(fc_target %>% dplyr::select(bin, FC, comparison), by = "bin")
    
    comparisons_lt[[i]] <- standard_lens_hex
  }
  
  # Combine all comparisons
  comparisons_df <- dplyr::bind_rows(comparisons_lt)
  
  return(comparisons_df)
}