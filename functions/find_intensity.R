find_intensity <- function (pixel_info = pixel_info, 
                            bin_data = bin_data, 
                            mz = mz, 
                            run = "all" ){
  library(stringr)
  
  df <- setNames(as.data.frame(str_split_fixed(colnames(bin_data), pattern = "@", n = 3)), c("run", "class", "bin"))
  df$bin <- as.numeric(df$bin)
  
  if (identical(run, "all")){
    df_mz <- bin_data[which(as.numeric(rownames(bin_data)) == as.numeric(mz)),]
    df$intensity <-  as.numeric(df_mz)
    pixel_intensity <- inner_join(pixel_info, df, by = c("run", "bin"))
    return(pixel_intensity)
  } else {
    index <- which(df$run %in% run) 
    df_run_flt <- df[index,]
    df_mz <- bin_data[which(as.numeric(rownames(bin_data)) == as.numeric(mz)), index]
    df_run_flt$intensity <- as.numeric(df_mz)
    pixel_intensity <- inner_join(pixel_info, df_run_flt, by = c("run", "bin"))
    return(pixel_intensity)
  }
}