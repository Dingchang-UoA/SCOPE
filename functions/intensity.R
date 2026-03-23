 intensity <- function (data = combinedimdata, 
                        mz, 
                        tol = tol, 
                        unit = unit){
   
   if (class(data)[1] == "MSImagingExperiment") {
     imdata = data
     ref <- mz(imdata)
     matched_mz <- match_mz(x = mz, ref = ref, tol = tol, unit = unit)
     matched_mz_idx <- which(mz(imdata) == matched_mz)
     
     imdata_flt <- imdata[matched_mz_idx,]
     sp <- Cardinal::spectra(imdata_flt, "intensity")
     mz_intensity <- as.numeric(sp[1, ]) 
     return(mz_intensity)
     
   } else if (class(data) == "data.frame") {
     
     suppressMessages(suppressWarnings(library(stringr)))
     bin_data <- data
     
     ref <- as.numeric(rownames(bin_data))
     matched_mz <- match_mz(x = mz, ref = ref, tol = tol, unit = unit)
     matched_mz_idx <- which(ref == matched_mz)
     
     df <- setNames(as.data.frame(str_split_fixed(colnames(bin_data), pattern = "@", n = 2)), c("run",  "bin"))
     df$bin <- as.numeric(df$bin)
     df_mz <- bin_data[matched_mz_idx,]
     df$intensity <-  as.numeric(df_mz)
     
     return(df)
     
   } else{
     stop("The data needs to be either imdata or bin data.")
   }
   
   
 }