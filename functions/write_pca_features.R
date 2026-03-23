write_pca_features <- function(PCA.contri, top_bin) {
  
  dimname <- colnames(PCA.contri)[!(colnames(PCA.contri) %in% c("mz", "class"))]
  
  Dim.feature.contri.list<- list()
  Dim.feature.contri.summmary.list<- list()
  
  for (dim in dimname){
    PCA.contri[,c("mz", "class",dim)]->dim_df
    colnames(dim_df)<-c("mz", "class","Contribution")
    dim_df$logcontrib<-log(dim_df$Contribution)
    dim_df$mz<-as.numeric(dim_df$mz)
    dim_df  %>% mutate(new_bin = ntile(Contribution, n=top_bin))->dim_df
    dim_df$Dim=dim
    write.csv(dim_df[dim_df$new_bin %in% c(top_bin),],paste0(dataset_path, "PCA_res/",dim,"_anno_freq_top1bin",".csv"),row.names = F)
    
    feature_slt<- unique(dim_df[dim_df$new_bin %in% c(top_bin),"mz"])
    Dim.feature.contri.list[[dim]]<-feature_slt
    Dim.feature.contri.summmary.list[[dim]]<-dim_df
  }
  Dim.feature.contri.summmary.list.bind<-do.call(rbind,Dim.feature.contri.summmary.list)
  
  Dim.feature.contri.summmary.list.bind %>% group_by(mz, Dim) %>% mutate(Rank = rank(Contribution)) -> Dim.feature.contri.summmary.list.bind
  
}


