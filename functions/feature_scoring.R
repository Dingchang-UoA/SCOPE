feature_scoring <- function (cpd_list = cpd_list) {
  
  library(HiTMaP)
  library(stringr)
  # setwd(paste0(function_path, "SIL_cpd_debug/"))
  # source(paste0(function_path, "SIL_cpd_debug/hitmap-wrapper.r"))
  # 
  raw_wd <- paste0(dataset_path) 
  fileinfo_v2 <- read.csv(paste0(dataset_path, "fileinfo_v2.csv"))
  locs_file<-fileinfo_v2$filenames
  # locs<-paste0(raw_wd,locs_file," ID/")
  cpd_list <- read.csv(paste0(dataset_path,"PCA_res/Metabo_Dim.1/mummichog_matched_compound_all.csv"))
  cmpd.map<-read.csv(paste0(dataset_path,"cmpd.map_anno.csv"))  
  cpd_list$formula<-cmpd.map$formula[match(cpd_list$Matched.Compound,cmpd.map$kegg_id)]
  cpd_list$Isotype<-NULL
  cpd_list<-cpd_list[!is.na(cpd_list$formula),]
  
  for (score_method in c("SQRTP")){
    for (timepoint in unique(locs_file)){
      #for (timepoint in "20hrSILGlucose_03_20211123"){
      #"16hrSILGlucose_01_20211123_trim"->timepoint
      for (FDR_cutoff in c(0.1))  {
        dir.create(paste0(dataset_path,score_method,"_",FDR_cutoff))
        peaklist<-as.data.frame(readr::read_csv(paste0(locs<-paste0("V:/Bioinformatics/imzML/",timepoint," ID/"),"overall_Spectrum.csv")))
        library(HiTMaP)  
        library(BiocParallel)
        HiTMaP:::Cpd_spectrum_match_rescore(cpd_list,peaklist,
                                            wd=paste0(dataset_path,"/",score_method,"_",FDR_cutoff),
                                            adjust_score=T,
                                            SPECTRUM_batch=timepoint,
                                            BPPARAM=SerialParam(),
                                            atom_isotope_sub=NULL,
                                            ppm=10,
                                            score_method=score_method,
                                            FDR_cutoff=FDR_cutoff,Decoy_search=T)
      }
    }
  }
  
  
  peak_res_sum<-NULL
  peak_res_sum_2nd<-NULL
  for (score_method in c("SQRTP")){
    dir.create(paste0(wd,"/",score_method))
    for (timepoint in unique(locs_file)){
      for (FDR_cutoff in c(0.1))  {
        dir.create(paste0(wd,"/",score_method,"_",FDR_cutoff))
        peak_res<-read.csv(paste0(wd,"/",score_method,"_",FDR_cutoff,"/",timepoint,"/","CPD_1st_ID_score_rank_SQRTP.csv"))
        peak_res2nd<-read.csv(paste0(wd,"/",score_method,"_",FDR_cutoff,"/",timepoint,"/","CPD_2nd_ID_score_rank_SQRTP.csv"))
        if(nrow(peak_res)>=1) peak_res$source<-paste0(score_method,"_",FDR_cutoff,"_",timepoint)
        peak_res_sum[[paste0(score_method,"_",FDR_cutoff,"_",timepoint)]]<-peak_res
        if(nrow(peak_res2nd)>=1) peak_res2nd$source<-paste0(score_method,"_",FDR_cutoff,"_",timepoint)
        peak_res_sum_2nd[[paste0(score_method,"_",FDR_cutoff,"_",timepoint)]]<-peak_res2nd
      }
    }
  }
  peak_res_sum_bind<-do.call(rbind,peak_res_sum)
  peak_res_sum_2nd_bind<-do.call(rbind,peak_res_sum_2nd)
  peak_res_sum_bind$CPD_file<-paste0(peak_res_sum_bind$formula,"_",
                                     peak_res_sum_bind$isdecoy,"_",
                                     peak_res_sum_bind$Isotype,"_",
                                     peak_res_sum_bind$source)
  peak_res_sum_2nd_bind$CPD_file<-paste0(peak_res_sum_2nd_bind$formula,"_",
                                         peak_res_sum_2nd_bind$isdecoy,"_",
                                         peak_res_sum_2nd_bind$Isotype,"_",
                                         peak_res_sum_2nd_bind$source)
  
  
  target_fm<-c("C6H12O6","C12H22O11","C15H24N2O17P2","C6H14O6","C7H15O10P","C6H13O9P1","C6H14O12P2")
  target_nm<-c("Glucose","Sucrose","UDP_Glucose","Sorbitol","S7P","G6P","F16P")
  peak_res_sum_bind_tgt<-peak_res_sum_bind[peak_res_sum_bind$formula_mono %in% target_fm,]
  peak_res_sum_bind_tgt_meta<-match(str_split_fixed(peak_res_sum_bind_tgt$source,"_",3)[,3],fileinfo_v2$filenames)
  peak_res_sum_bind_tgt<-cbind(peak_res_sum_bind_tgt,fileinfo_v2[peak_res_sum_bind_tgt_meta,])
  peak_res_sum_bind_tgt$FDR_pass<-peak_res_sum_bind_tgt$CPD_file %in% unique(peak_res_sum_2nd_bind$CPD_file)
  peak_res_sum_bind_tgt<-unique(peak_res_sum_bind_tgt[,c("Score" ,"adduct", "formula","formula_mono","Isotype","isdecoy","mins","bioreplicate","class" ,"FDR_pass" )])
  peak_res_sum_bind_tgt$cpd<-target_nm[match(peak_res_sum_bind_tgt$formula_mono,target_fm)]
  
  data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sem = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
  }
  library(ggplot2)
  df2 <- data_summary(peak_res_sum_bind_tgt, varname="Score", 
                      groupnames=c("cpd","Isotype","isdecoy","class","adduct"))
  df2 <- merge(df2,peak_res_sum_bind_tgt[,c("cpd","Isotype","isdecoy","class","adduct","FDR_pass")],by=c("cpd","Isotype","isdecoy","class","adduct"))
  p<-ggplot(df2) +geom_point(aes(x=class,y=Score,colour=FDR_pass,size=FDR_pass)) + geom_line(aes(x=class,y=Score,color=interaction(isdecoy))) + 
    geom_errorbar(aes(x=class,ymin=Score-sem, ymax=Score+sem), width=.2,
                  position=position_dodge(.9)) + 
    facet_grid(rows = vars(Isotype,adduct),cols = vars(cpd))
  png(paste0(wd,"/","/target_fm_summary_",FDR_cutoff,".png"),height = 2400*2,width = 1500*length(target_nm),res = 300)
  
  print(p)
  dev.off()    
  write.csv(peak_res_sum_bind,paste0(wd,"/",score_method,"_","peak_res_sum_bind_",FDR_cutoff,".csv"),row.names = F)
  write.csv(peak_res_sum_2nd_bind,paste0(wd,"/",score_method,"_","peak_res_sum_2nd_bind_",FDR_cutoff,".csv"),row.names = F)
}

