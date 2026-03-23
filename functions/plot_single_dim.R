plot_single_dim <- function(dims, group, raw_info = pixel_info, run = pixel_info$run, fig_path){
   
  dims <- as.data.frame(dims)
  colnames(raw_info)[3] <- c("group")
  
  for (i in 1: ncol(dims) ){
    Rank <- dense_rank(dims[,i])
    df <- as.data.frame(cbind(dims[,i], Rank, group))
    
    raw_info$Rank <- as.numeric(df$Rank[match(raw_info$group, df$group)])
    raw_info$run <- run
    
    colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(max(raw_info$Rank))
    
    png(paste0(fig_path, "tsne_dim_2.png"), units = "px", width = 9000, height = 6000, res = 300)
    ggplot(raw_info, aes(z.final, r.final, color = Rank))+
      geom_point()+
      scale_color_gradientn(colors = colors) +  
      facet_wrap(vars(run))+ 
      egg::theme_article()+
      coord_equal()
    dev.off()
  }
}