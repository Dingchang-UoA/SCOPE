virtual_segmentation <- function(coordatadf,
                                Virtual_segmentation_rankfile = Virtual_segmentation_rankfile){
  library(reshape2)
  library(BiocParallel)
  library(data.table)
  
  coordatadf <- as.data.frame(coordatadf)
  
  coordatadf$indx<-as.numeric(rownames(coordatadf))
  coordatadf<-as.data.frame(coordatadf)
  
  coordatalist=list()
  
  for (i in unique(coordatadf$run)){
    
    coordatalist[[i]]<-coordatadf[coordatadf$run==i,]
    
  }
  
  
  radius_rank=read.csv(file = Virtual_segmentation_rankfile)
  radius_rank=radius_rank[order(radius_rank$Rank),]
  
  
  coordata_rank_list<- function(coordata,radius_rank,BPPARAM=bpparam()){
    coordist_para = function(i,coordata){
      
      coordist_para=NULL
      for( j in 1  :  nrow(coordata)){
        coordist_para[j]=sqrt((coordata$x[i]-coordata$x[j])^2+(coordata$y[i]-coordata$y[j])^2)
      }
      coordist_para
    }
    
    
    #coordistmatrix<-matrix(nrow=nrow(coordata),ncol=nrow(coordata))
    #coordistmatrix<-NULL
    #cl=autoStopCluster(cl)
    #coordistmatrix=parLapply(cl=cl,1:  nrow(coordata),coordist_para,coordata)
    coordistmatrix=bplapply(1:  nrow(coordata),coordist_para,coordata,BPPARAM = BPPARAM)
    coordistmatrix=matrix(unlist(coordistmatrix),nrow = nrow(coordata),ncol = nrow(coordata))
    coordistmatrix=as.data.table(coordistmatrix)
    coordistmatrix$sum=0
    #coordistmatrix$sum=base::unlist(parLapply(cl=cl,1:nrow(coordata),function(j,coordistmatrix,coordata){coordistmatrix$sum[j]=sum(coordistmatrix[j,1:nrow(coordata)])},coordistmatrix,coordata))
    coordistmatrix$sum=base::unlist(bplapply(1:nrow(coordata),function(j,coordistmatrix,coordata){coordistmatrix$sum[j]=sum(coordistmatrix[j,1:nrow(coordata)])},coordistmatrix,coordata,BPPARAM = BPPARAM))
    
    coorrange=max(coordistmatrix$sum)-min(coordistmatrix$sum)
    
    
    
    findedge<-function(coordata){
      #for (i in 1: nrow(coordata)){
      uniquex=unique(coordata$x)
      uniquey=unique(coordata$y)
      coordata$edge=FALSE
      for (x in uniquex){
        
        min=min(coordata[coordata$x==x,"y"])
        max=max(coordata[coordata$x==x,"y"])
        coordata['&'(coordata$y==max,coordata$x==x),"edge"]=TRUE
        coordata['&'(coordata$y==min,coordata$x==x),"edge"]=TRUE
      }
      
      for (y in uniquey){
        min=min(coordata[coordata$y==y,"x"])
        max=max(coordata[coordata$y==y,"x"])
        coordata['&'(coordata$x==max,coordata$y==y),"edge"]=TRUE
        coordata['&'(coordata$x==min,coordata$y==y),"edge"]=TRUE
      }
      coordata
    }
    
    coordata=findedge(coordata)
    
    rank_pixel<-function(coordata,coordistmatrix){
      shape_center=coordata[coordistmatrix$sum==min(coordistmatrix$sum),]
      center_dist=t(coordistmatrix[which.min(coordistmatrix$sum),1:nrow(coordata)])
      library(useful)
      From <- shape_center[rep(seq_len(nrow(shape_center)), each=nrow(coordata)),1:2]
      To <- coordata[,1:2]
      df=To-From
      center_edge_angle<-cbind(coordata[,1:2],cart2pol(df$x, df$y, degrees = F),edge=coordata[,"edge"])
      center_edge_angle_sdge=center_edge_angle[center_edge_angle$edge==TRUE,]
      coordata$rank=0
      coordata$pattern=""
      
      for (i in 1: (nrow(coordata))){
        
        if (coordata$edge[i]!=TRUE){
          df=coordata[i,1:2]-shape_center[,1:2]
          point_center_angle<-cbind(coordata[i,1:2],cart2pol(df$x, df$y, degrees = F))
          pointedge=center_edge_angle_sdge[which(abs(center_edge_angle_sdge$theta-point_center_angle$theta)==min(abs(center_edge_angle_sdge$theta-point_center_angle$theta))),]
          
          pointedge=pointedge[which.min(pointedge$r),]
          to_edge=coordistmatrix[[i]]['&'(coordata$x==pointedge$x,coordata$y==pointedge$y)]
        }else{to_edge=0}
        
        
        to_center=center_dist[i]
        total=to_edge+to_center
        max(radius_rank$Radius_U)
        norm_center_dist=to_center/total*max(radius_rank$Radius_U)
        
        
        if (norm_center_dist == 0) {
          norm_center_dist <- norm_center_dist + 0.00000000000000000000000001
          coordata$rank[i]=as.character(radius_rank$Rank['&'(radius_rank$Radius_L<norm_center_dist,radius_rank$Radius_U>=norm_center_dist)])
          coordata$pattern[i]=as.character(radius_rank$Name['&'(radius_rank$Radius_L<norm_center_dist,radius_rank$Radius_U>=norm_center_dist)])
          
        } 
        
        coordata$rank[i]=as.character(radius_rank$Rank['&'(radius_rank$Radius_L<norm_center_dist,radius_rank$Radius_U>=norm_center_dist)])
        coordata$pattern[i]=as.character(radius_rank$Name['&'(radius_rank$Radius_L<norm_center_dist,radius_rank$Radius_U>=norm_center_dist)])
      }
      coordata
    }
    
    coordata=rank_pixel(coordata,coordistmatrix)
    
    coordata
    
  }
  
  
  coordatalist_result<-lapply(coordatalist,coordata_rank_list,radius_rank, BPPARAM = SnowParam(workers = 20, progressbar = T))
  
  coordatalist_result_merge<-do.call(rbind,coordatalist_result)
  coordatadf_merge<-merge(coordatadf,coordatalist_result_merge[,c("x","y","run","rank","pattern")],by=c("x","y","run"))
  rownames(coordatadf_merge)<-coordatadf_merge$indx
  return(coordatadf_merge)
  
}
