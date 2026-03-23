Plot_FC_HeatMap <- function(hexbin_average_df = hexbin_average_df,
                            comparisons = comparisons,
                            mz = mz, 
                            HexBin_index = HexBin_index,
                            groups = c("hbnc", "hbnn", "hboc", "hbon")
                              ){
  
  comparisons_lt <- list()
  
  for (comparison in 1:length(comparisons) ) {
    hexbin_average_df <- hexbin_average_df
    info <- as.data.frame(str_split_fixed(colnames(hexbin_average_df), "@", 2))
    colnames(info) <- c("treatment", "bin")
    
    hexbin_average_flt <- hexbin_average_df[, which(info$treatment %in% comparisons[[comparison]])]
    
    hbn_mean_raw <- hexbin_average_flt %>% select(matches(paste0("^", comparisons[[comparison]][1])))
    hbo_mean_raw <- hexbin_average_flt %>% select(matches(paste0("^", comparisons[[comparison]][2])))
    hbn_mean_raw <- hbn_mean_raw+1
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
    fc_bin_raw_mat <- as.matrix(fc_bin_raw_df)
    
    na_idx <- as.data.frame(which(is.na(fc_bin_raw_mat), arr.ind = TRUE))
    
    rownames(fc_bin_raw) <- fc_bin_raw$mz
    fc_bin_df <- fc_bin_raw[,-1]
    fc_bin_mat <- as.matrix(fc_bin_df)
    
    # 
    fc_bin_mat[as.matrix(na_idx)] <- fc_bin_raw_mat[as.matrix(na_idx)]
    
    fc_bin_clean <- as.data.frame(fc_bin_mat)
    rownames(fc_bin_clean) <- rownames(fc_bin_raw)
    
    
    fc_GSH <- as.data.frame(t(fc_bin_clean[which(rownames(fc_bin_clean) == mz),]))
    colnames(fc_GSH) <- c("FC")
    
    
    fc_GSH$bin <- as.numeric(rownames(fc_GSH))
    fc_GSH$comparison <- paste0(comparisons[[comparison]][1], " vs " , comparisons[[comparison]][2] )
    
    standard_lens_hex <- standard_lens %>% 
      left_join(fc_GSH %>% select(bin, FC, comparison),
                by = c( "bin"))
    
    comparisons_lt[[comparison]] <- standard_lens_hex
    
    
  }
  
  comparisons_df <- do.call(rbind, comparisons_lt)
  
  
  bin_916 <- comparisons_df %>%
    filter(bin == HexBin_index)
  
  df_heat <- bin_916  %>%
    select(comparison, FC)
  
  df_heat <- df_heat %>%
    separate(comparison, into = c("group1", "group2"), sep = " vs ")
  
  
  df_heat <- bin_916 %>%
    tidyr::separate(comparison, into = c("group1", "group2"), sep = " vs ")
  
  df_sym <- bind_rows(
    df_heat,
    df_heat %>% transmute(group1 = group2, group2 = group1, FC = FC)
  )
  
  
  df_diag <- data.frame(
    group1 = groups,
    group2 = groups,
    FC = 1
  )
  
  df_all <- bind_rows(df_sym, df_diag)
  
  df_all <- df_all %>%
    mutate(group1 = factor(group1, levels = groups),
           group2 = factor(group2, levels = groups))
  
  p <- ggplot(df_all, aes(x = group1, y = group2, fill = log2(FC))) +
    geom_tile(color = "white") +
    scale_fill_distiller(
      palette = "RdYlBu",
      type = "div",
      direction = -1,
      limits = c(-3, 3),
      oob = scales::squish,
      name = expression(log[2](FC))
    ) +
    coord_fixed() +
    theme_article() +
    theme(legend.position = "inside", 
          legend.position.inside  = c(0.88, 0.2),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  return(p)
}