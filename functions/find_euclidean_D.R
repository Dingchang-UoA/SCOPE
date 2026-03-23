find_euclidean_D <- function(bin_data){
  
  library(stringr)
  library(stats)
  library(dplyr)
  
  # Extract metadata from column names
  df <- setNames(as.data.frame(str_split_fixed(colnames(bin_data), pattern = "@", n = 3)), 
                 c("run", "class", "bin"))
  df$bin <- as.numeric(df$bin)
  df$class <- as.numeric(df$class)
  
  euclidean_list <- list()
  
  for (class in unique(df$class)){
    df_flt <- df[which(df$class == class),] 
    unique_bins <- unique(df_flt$bin)
    unique_runs <- unique(df_flt$run)
    # Initialize results matrix with an extra column for bin index
    euclidean_results <- matrix(NA, 
                                nrow = length(unique_bins), 
                                ncol = 4)
    colnames(euclidean_results) <- c("bin", 
                                     "Sample1_vs_Sample2",
                                     "Sample1_vs_Sample3",
                                     "Sample2_vs_Sample3")
                                    
    euclidean_results[, 1] <- unique_bins 
    
    
    for (bin_idx in unique_bins){
      index <- which(df$class == class & df$bin == bin_idx)  
      index <- index[order(match(df$class[index], class))] 
      bin_data_flt <- bin_data[,index, drop = FALSE]
      mz_bin_df <- matrix(0, nrow = nrow(bin_data_flt), ncol = 3)
      colnames(mz_bin_df) <- unique_runs 
      
      col.nms.split <- as.data.frame(str_split_fixed(colnames(bin_data_flt), "@", 3))
      colnames(col.nms.split) <- c("run", "class", "bin")
      
      for ( i in col.nms.split$run){
        mz_bin_df[,i] <- bin_data_flt[,which(col.nms.split$run == i)]
      }
      
      sample1 <- as.vector(mz_bin_df[, 1, drop = FALSE])
      sample2 <- as.vector(mz_bin_df[, 2, drop = FALSE])
      sample3 <- as.vector(mz_bin_df[, 3, drop = FALSE])
      
      e1 <- as.numeric(dist(rbind(sample1, sample2)))
      e2 <- as.numeric(dist(rbind(sample1, sample3)))
      e3 <- as.numeric(dist(rbind(sample2, sample3)))
      
      e_df <- data.frame(e1,e2,e3)
      colnames(e_df) <- c("Sample1_vs_Sample2",
                          "Sample1_vs_Sample3",
                          "Sample2_vs_Sample3")
      rownames(e_df) <- bin_idx
      
      for (sample in colnames(e_df)){
        euclidean_results[which(euclidean_results[, c("bin")] == as.numeric(bin_idx)), sample] <-  e_df[as.character(bin_idx),sample]           
      }
      
      euclidean_results <- as.data.frame(euclidean_results)
      euclidean_results <- euclidean_results %>%
        mutate(dist_mean = rowMeans(select(., Sample1_vs_Sample2, Sample1_vs_Sample3, Sample2_vs_Sample3), na.rm = TRUE),
               class = class)
      
    }
    
    euclidean_list[[as.character(class)]] <- euclidean_results
  }
  
  df_dist <- do.call(rbind, euclidean_list)

  return(df_dist)  # Convert to dataframe for better readability
}
