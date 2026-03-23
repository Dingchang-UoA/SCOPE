plot_pip_fig <- function ( raw_data = cls_merge_bind, method = c("UMAP", "t-SNE", "PCA"), group = pixel_info$label){
  
  set.seed(2025)
  library(stringr)
  col.nms.split<-as.data.frame(str_split_fixed(rownames(raw_data),"@",2))
  colnames(col.nms.split) <- c("mz", "class")
  
  umap_list <-list()
  
  for (class in unique(col.nms.split$class)){
    raw_data_class <- raw_data[ which(col.nms.split$class == class) ,]
    raw_data_class_mt <- t(raw_data_class)
    raw_data_class_mt[is.na(raw_data_class_mt)] <- 0
    set.seed(42)
    features.umap  <-uwot::umap(raw_data_class_mt,
                                n_neighbors = 15,
                                n_components = 2, 
                                metric = "euclidean",
                                min_dist = 0.1,
                                spread = 1,
                                local_connectivity = 1,
                                repulsion_strength = 1,
                                negative_sample_rate = 5,
                                learning_rate = 1)
    
    
    features.umap <- as.data.frame(features.umap)
    colnames(features.umap) <- c("x", "y")
    features.umap$bin <- rownames(features.umap)
    features.umap$class <- class
    umap_list[[as.character(class)]] <- features.umap
  }
  umap_df <- do.call(rbind, umap_list)
  umap_df$label <- pixel_info$label[match(as.numeric(umap_df$bin), as.numeric(pixel_info$bin))]
  
}