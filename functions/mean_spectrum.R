mean_spectrum <- function(imdata, run = c("all")){
  
  if (run == "all") {
    
    feature_intensity <- as.data.frame(as.matrix(spectra(imdata, "intensity")))
    rownames(feature_intensity) <- mz(imdata) 
    feature_intensity_mean <- data.frame(mz = rownames(feature_intensity), 
                                         mean_spectrum = rowMeans(feature_intensity))
    
  } else{
    imdata_flt <- combinedimdata_raw[, which(as.character(run(imdata)) == run )]
    
    feature_intensity <- as.data.frame(spectra(imdata_flt, "intensity"))
    rownames(feature_intensity) <- mz(imdata) 
    feature_intensity_mean <- data.frame(mz = rownames(feature_intensity), 
                                         mean_spectrum = rowMeans(feature_intensity))
  }
  
  return(feature_intensity_mean)
}