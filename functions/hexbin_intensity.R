hexbin_intensity <- function(imdata, pixel_info ){
  
  combinedimdata<-imdata
  
  suppressMessages(suppressWarnings(library(hexbin)))
  
  label_sel="all"
  merged_data_df<-NULL
  coordata_label_sel_arc<-NULL
  
  
  setCardinalBPPARAM(BPPARAM=SerialParam())
  setCardinalVerbose(verbose=F)
  
  for (imsrun in unique(pixel_info$run)){
    message(imsrun)
    coordata_run<-pixel_info[pixel_info$run==imsrun,]
    rownames(coordata_run)<-1:nrow(coordata_run)
    
    combinedimdata_sample_bin <- combinedimdata[,run(combinedimdata)==imsrun]
    
    BiocParallel::bplapply(unique(coordata_run$bin),function(bincell,coordata_run,combinedimdata_sample_bin){
      meanspec <- combinedimdata_sample_bin [,which(coordata_run$bin==bincell)] %>% 
        summarizeFeatures( FUN = "mean") %>% 
        process()
      return(meanspec@featureData@listData[["mean"]])
    },coordata_run,combinedimdata_sample_bin,BPPARAM = SnowParam(workers = 2))->temp
    
    
    
    bincell_df_bind<-do.call(cbind,temp)
    colnames(bincell_df_bind)<-paste(imsrun, as.character(unique(coordata_run$bin)),sep = "@")
    merged_data_df[[imsrun]]<-bincell_df_bind
    coordata_label_sel_arc[[imsrun]]<-coordata_run
  }
  
  coordata_label_sel_arc_bind<-do.call(rbind,coordata_label_sel_arc)
  merged_data_df_bind<-do.call(cbind,merged_data_df)
  rownames(merged_data_df_bind)<-combinedimdata@featureData@mz
  
  #Generate hexbin data
  write.csv(merged_data_df_bind,paste0(dataset_path,"/",label_sel,"_merged_data_df_bind",".csv"))
  
  #Generate Master Table 
  write.csv(coordata_label_sel_arc_bind,paste0(dataset_path,"/",label_sel,"_coordata_label_sel_arc_bind",".csv"))
  
  
  message(label_sel," DONE")
}
