MSI_preprocessing <- function(datafile, workdir, deconv_peaklist, ppm, import_ppm = ppm/2, threshold = 1e-07, rotate = NULL) {
  
  suppressMessages(suppressWarnings(library(matter)))
  suppressMessages(suppressWarnings(library(stringr)))
  suppressMessages(suppressWarnings(library(HiTMaP)))
  suppressMessages(suppressWarnings(library(Cardinal)))
  datafile_base <- basename(datafile)
  datafile <- str_remove(datafile_base, "\\.imzML$")
  
  if (length(workdir) != length(datafile)) {
    workdir = rep(workdir[1], length(datafile))
  }
  
  datafile_imzML = paste0(datafile, ".imzML")
  
  deconv_peaklist_decov <- HiTMaP:::isopattern_ppm_filter_peaklist(deconv_peaklist, 
                                                                   ppm = ppm, threshold = threshold)
  
  imdata <<- NULL
  imdata <- NULL
  gc()
  combinedimdata_list <- NULL
  
  
  for (z in 1:length(datafile)) {
    deconv_peaklist_decov_plot <- deconv_peaklist_decov
    if (!("m.z" %in% colnames(deconv_peaklist_decov_plot))) {
      deconv_peaklist_decov_plot$m.z <- deconv_peaklist_decov_plot[, 1]
    }
      
    setwd(workdir[1])
    message("Loading raw image data for preprocessing: ", 
              paste0(gsub(".imzML$", "", datafile[z]), ".imzML"))
      
    imdata <- Cardinal::readMSIData(paste0(workdir[z], 
                                               "/", datafile_imzML[z]), resolution = import_ppm, 
                                        units = "ppm", BPPARAM = SerialParam(), memory = T) %>%  process()
           
     
    if (!is.null(rotate)){
      
      if (rotate[z] == 0){
        message("No need for rotation.")
      } else {
        imdata <- rotateMSI_2D(imdata = imdata, rotation_degree = rotate[z])
        message(paste0("imdata has been rotated for ", rotate[z], " degree."))
      }
    }
    
    
   
    message("mzbin for ", datafile[z])
    setCardinalBPPARAM(SerialParam())
    
    New_fdata_LB <- imdata[1, ]
    New_fdata_LB@featureData@listData[["mz"]][1] <- min(deconv_peaklist_decov_plot$m.z, 
                                                        na.rm = T) - 2 * ppm/1e+06 * min(deconv_peaklist_decov_plot$m.z, 
                                                                                         na.rm = T)
    New_fdata_UB <- New_fdata_LB
    New_fdata_UB@featureData@listData[["mz"]][1] <- max(deconv_peaklist_decov_plot$m.z, 
                                                        na.rm = T) + 2 * ppm/1e+06 * min(deconv_peaklist_decov_plot$m.z, 
                                                                                         na.rm = T)
    if (min(mz(imdata)) > min(New_fdata_LB@featureData@listData[["mz"]][1])) {
      message("removing ", sum(deconv_peaklist_decov_plot$m.z < 
                                 min(mz(imdata))), " m/z features from target list, out of LB of ", 
              datafile[z])
      deconv_peaklist_decov_plot <- deconv_peaklist_decov_plot[deconv_peaklist_decov_plot$m.z >= 
                                                                 min(mz(imdata)), ]
      deconv_peaklist_decov <- deconv_peaklist_decov_plot
    }
    if (max(mz(imdata)) < max(New_fdata_UB@featureData@listData[["mz"]][1])) {
      message("removing ", sum(deconv_peaklist_decov_plot$m.z > 
                                 max(mz(imdata))), " m/z features from target list, out of UB of ", 
              datafile[z])
      deconv_peaklist_decov_plot <- deconv_peaklist_decov_plot[deconv_peaklist_decov_plot$m.z <= 
                                                                 max(mz(imdata)), ]
      deconv_peaklist_decov <- deconv_peaklist_decov_plot
    }
    imdata_bin <- imdata %>% bin(ref =unlist(deconv_peaklist_decov_plot$m.z), 
                                 method = "sum",
                                   tolerance = ppm/2,
                               units = "ppm") %>% process()
    
  
    #imdata@elementMetadata@coord@listData[["z"]] <- NULL
    #imdata@elementMetadata@resolution = c(x = 1, y = 1)
    combinedimdata_list[[z]] <- imdata_bin
    rm(imdata)
    rm(imdata_bin)
  }
  
  for (z in 1:length(datafile)) {
    combinedimdata_list[[z]] <- combinedimdata_list[[z]][mz(combinedimdata_list[[z]]) %in% 
                                                           deconv_peaklist_decov_plot$m.z, ]
  }
  gc()
  saveRDS(combinedimdata_list, paste0(workdir[1], "/", "combinedimdata_list.rds"), compress = T)
  for (z in 1:length(datafile)) {
    combinedimdata_list[[z]]@centroided <- TRUE
    if (z == 1) {
      combinedimdata <- combinedimdata_list[[z]]
    }
    else {
      combinedimdata <- Cardinal::cbind(combinedimdata, 
                                        combinedimdata_list[[z]])
    }
    combinedimdata_list[[z]] <- 1
    gc()
  }
  saveRDS(combinedimdata, paste0(workdir[1], "/combinedimdata.rds"), compress = T)
  return(paste0(workdir[1], 
                "/combinedimdata.rds"))
  
}