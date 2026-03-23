compare_groups_stats <-  function (mean_sp_df, mz, comparisons, regions = c("OC", "IC", "Core") ){
  info <- as.data.frame(stringr::str_split_fixed(colnames(mean_sp_df), "@", 2))
  colnames(info) <- c("run", "region")
  info$treatment = sub("_.*", "", info$run)
  mz_idx <- which(rownames(mean_sp_df) == mz)
  
  result_list <- list()
  
  for (comparison in 1: length(comparisons)) {
    group1 <- comparisons[[comparison]][1]
    group2 <- comparisons[[comparison]][2]
    
    for (region in regions){
      idx_group1 <- info$treatment == group1
      idx_group2 <- info$treatment == group2
      idx_region <- info$region == region
      
      group1_values <- unlist(mean_sp_df[mz_idx, '&'(idx_group1, idx_region)])
      group2_values <- unlist(mean_sp_df[mz_idx, '&'(idx_group2, idx_region)])
      
      if (var(group1_values, na.rm = TRUE) == 0 && var(group2_values, na.rm = TRUE) == 0) {
        p_value <- NA
        log2fc <- NA
      } else {
        t_res <- t.test(group1_values[!is.na(group1_values)],
                        group2_values[!is.na(group2_values)])
        p_value <- t_res$p.value
        log2fc <- log2(mean(group2_values, na.rm = TRUE) / mean(group1_values, na.rm = TRUE))
      }
      result_list[[length(result_list) + 1]] <- data.frame(
        m.z = mz,
        group1 = group1,
        group2 = group2,
        region = region,
        log2FC = log2fc,
        p_value = p_value
      )
      
      
    }
    
  }
  return(do.call(rbind, result_list))
}