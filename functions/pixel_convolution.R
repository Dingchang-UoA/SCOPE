pixel_convolution <- function(umap.res, 
                              sp,
                              off_tissue_cluster, 
                              convolution_level = c("square_bin", "region")[1]){
  
  off_tissue_idx <- which(umap.res$cluster %in% off_tissue_cluster)
  sp_on_tissue <- sp
  sp_on_tissue[, off_tissue_idx] <- NA 
    
  if (convolution_level == "region"){
    
    mean_sp_list <- list()
    
    for (region in unique(umap.res$pattern)) {
      idx <- which(umap.res$pattern == region)
      sp_flt <- sp_on_tissue[, idx]
      sp_mean <- rowMeans(sp_flt, na.rm = TRUE)
      
      mean_sp_list[[paste0(unique(umap.res$run), "@", region)]] <- sp_mean 
    }
    
    mean_sp_df <- do.call(cbind, mean_sp_list)
    mean_sp_df <- as.data.frame(mean_sp_df)
    
    return(mean_sp_df)
    
  } else if (convolution_level == "square_bin") {
    
    coordata_conv <- unique(umap.res[,c("x_conv", "y_conv", "pattern")])
    conv_sp_list <- list()
    
    for (i in 1:nrow(coordata_conv)){
      
      condition1 <- umap.res$x_conv == coordata_conv$x_conv[i]
      condition2 <- umap.res$y_conv == coordata_conv$y_conv[i]
      
      idx <- '&'(condition1, condition2)
      
      sp_flt <- as.data.frame(sp_on_tissue[,idx])
      sp_mean <- rowMeans(sp_flt, na.rm = T)
      
      conv_sp_list[[paste(unique(umap.res$run), 
                          coordata_conv$x_conv[i],
                          coordata_conv$y_conv[i], 
                          coordata_conv$pattern[i], sep = "@")]] <- sp_mean
      
    }
    
    mean_sp_df <- do.call(cbind, conv_sp_list)
    mean_sp_df <- as.data.frame(mean_sp_df)
    return(mean_sp_df)
    
  }
  
  
}