group_comparison <- function (data = bin_data_mean, pixel_info, condition, k_means_cluster){
  
  #Reference group goes first then experiment group
  
  pixel_group1 <- pixel_info %>% 
    filter(treatment == condition[1], cluster == k_means_cluster[1])
  
  pixel_group2 <- pixel_info %>% 
    filter(treatment == condition[2], cluster == k_means_cluster[2])
  
  pixel_group1$match <- paste(pixel_group1$treatment, pixel_group1$bin, sep = "@")
  
  pixel_group2$match <- paste(pixel_group2$treatment, pixel_group2$bin, sep = "@")
  
  idx_1 <- which(colnames(data) %in% pixel_group1$match)
  idx_2 <- which(colnames(data) %in% pixel_group2$match) 
  group_1 <- bin_data_mean[,idx_1]
  group_2 <- bin_data_mean[,idx_2]
  
  results <- list()
  
  for (mz in rownames(group_1)){
    values_1 <- as.numeric(unlist(group_1[mz, ]))
    values_2 <- as.numeric(unlist(group_2[mz, ]))
    t_res <- t.test(values_1, values_2)
    fc <- mean(values_2, na.rm = TRUE) / mean(values_1, na.rm = TRUE)
    
    results[[mz]] <- c(
      p.value = t_res$p.value,
      t.score = unname(t_res$statistic),
      FC = fc
    )
    
  }
  
  
  results <- do.call(rbind, results)
  results <- as.data.frame(results)
  results$m.z <- rownames(results)
  
  message(paste0("Remove ", sum(is.nan(results$p.value)), " features that have 0 variance.") )
  results <- results[!is.nan(results$p.value), ]
  
  return(results)
}