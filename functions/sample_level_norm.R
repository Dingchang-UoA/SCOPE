sample_level_norm <- function(imdata = combinedimdata, mz_ref, tol,
                              unit) {
  coord <- retrieve_coord(imdata)
  coord$intensity <- intensity(data = imdata, 
                               mz = mz_ref,
                               tol = tol,
                               unit = unit)

  library(dplyr)
  intensity_summary <- coord %>% group_by(run) %>%
    dplyr::summarise(mean_intensity = mean(intensity))
  
  
  run_list<-list()
  
  for (run in unique(run(imdata))){
    combinedimdata_subset<- imdata[, run(imdata) == run]
    message(unique(run(combinedimdata_subset)))
    intensity_summary$mean_intensity[which(intensity_summary$run == run)]->spl_mean
    combinedimdata_subset_norm <- addProcessing(combinedimdata_subset, FUN = function(x,mz){x / spl_mean} , label = "sample_norm")
    combinedimdata_subset_norm <- process(combinedimdata_subset_norm)
    run_list[[run]]<- combinedimdata_subset_norm
  }
  
  combinedimdata_norm_2 <-  do.call(cbind, run_list)
  
  return(combinedimdata_norm_2)
}



