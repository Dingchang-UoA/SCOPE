pixel_level_norm <- function(imdata = combinedimdata, mz_ref, tol,
                          unit) {
  
  coord <- NULL
  
  coord <- as.data.frame(imdata@elementMetadata@coord@listData)
  coord$run <- run(imdata)
  
  for (mz in mz_ref) {
    matched_mz <- match_mz(mz, mz(imdata), tol = tol, unit = unit)
    coord[, paste0(matched_mz)] <- intensity(data = imdata, 
                                 mz = mz,
                                 tol = tol,
                                 unit = unit)
      }
  
  matched_mz_all <- colnames(coord)[4:ncol(coord)]
  
  library(dplyr)
  intensity_summary <- coord %>% 
    select(-x, -y) %>%           
    group_by(run) %>%           
    summarise(across(everything(), mean, na.rm = TRUE)) %>% 
    mutate(mean_intensity = rowMeans(select(., -run)))  
  
    
  
  
  run_list<-list()
  
  for (run in unique(run(imdata))){
    combinedimdata_subset<- imdata[, run(imdata) == run]
    message(unique(run(combinedimdata_subset)))
    combinedimdata_subset_norm <- process(combinedimdata_subset, 
                                          function(x) x / intensity_summary$mean_intensity[which(intensity_summary$run == run)], 
                                          label = "Normalization")
    run_list[[run]]<- combinedimdata_subset_norm
  }
  
  combinedimdata_norm_2 <-  do.call(cbind, run_list)
  
  return(combinedimdata_norm_2)
}