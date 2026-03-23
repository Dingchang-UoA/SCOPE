write_pathway_enrichment_result <- function(Dims = paste0("Dim.", 1:5), outwd = outwd ) {
  
    for  (dim in Dims) {
    
    dim.sig<-readr::read_csv(paste0(outwd, dim, "_anno_freq_top1bin.csv"))
    dir.create(paste0(outwd,"/Metabo_",dim))
    setwd(paste0(outwd,"/Metabo_",dim))
    message(getwd())
    
    mSet<-NULL
    
    mSet<-InitDataObjects("mass_table", "mummichog", FALSE)
    anal.type<<-mSet$analSet$type
    mSet<-SetPeakFormat(mSet, "pvalue")
    mSet<-UpdateInstrumentParameters(mSet, 10.0, "negative", "yes");
    mSet<-SetRTincluded(mSet, "no")
    mSet<-Read.TextData(mSet, paste0(dataset_path, "/bin_data_time_sumspec_tb.csv"), "colu", "disc")
    mSet<-SanityCheckMummichogData(mSet)
    mSet<-ReplaceMin(mSet);
    mSet<-SanityCheckMummichogData(mSet)
    mSet<-PreparePrenormData(mSet)
    mSet<-Normalization(mSet, "none", "none", "none", ratio=FALSE, ratioNum=20)
    
    
    
    mSet<-SetPeakEnrichMethod(mSet, "mum", "v2")
    mSet<-PreparePeakTable4PSEA(mSet)
    
    feature_sig <-  names(mSet[["analSet"]][["tt"]][["t.score"]]) %in% dim.sig$mz
    mSet[["analSet"]][["tt"]][["p.value"]][feature_sig] <-0.01
    mSet[["analSet"]][["tt"]][["p.value"]][!feature_sig] <-1
    nData<-mSet$dataSet$mummi.proc
    feature_sig <-  nData$m.z %in% dim.sig$mz
    match_mz(nData$m.z,dim.sig$mz,10,"ppm")->mapped_mz
    feature_sig <- !is.na(mapped_mz)
    nData[feature_sig,1] <-0.000001
    nData[!feature_sig,1] <-1
    nData->mSet$dataSet$mummi.proc
    
    mSet<-SetMummichogPval(mSet, 0.05)
    mummichog.lib.new<-mummichog.lib.mod(lib="bta_kegg",lib.new="bta_kegg_13C_adduct",C13_number=c(6),C13_deltamass=1.003355,wd=paste0(outwd,"/Metabo_",dim),force_update_lib=F,adducts_list=c("M+Cl[-]","M-H[-]"),method = "new_adduct")  
    mSet$paramSet$ContainsMS2<-F
    mSet<-PerformPSEA(mSet, "bta_kegg_13C_adduct", "current", 3 , 100)
    mSet<-SaveTransformedData(mSet)
    mSet<-PlotPeaks2Paths(mSet, "peaks_to_paths_0_", "png", 300, width=NA)
    write.csv(annotate_PSEA(mSet = mSet,mummichog.lib.new),"mSet_mummi_resmat_anno.csv")
    
  }
  
  df_C13<-NULL
  for (p_thres in paste0("Metabo_Dim.",1:5)){
    df_C13[[as.character(p_thres)]]<-read.csv(paste0(outwd,"/",p_thres,"/","/","mSet_mummi_resmat_anno.csv"))
    df_C13[[as.character(p_thres)]]$source<-p_thres
  }
}


