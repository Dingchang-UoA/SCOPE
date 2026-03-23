Correlation_plot <- function(HexBin_FC = HexBin_FC,
                            pixel_info = pixel_info,
                            mz1, 
                            mz2,
                            umap.res.ref = umap.res.ref){
  
  library(patchwork)
  library(ggpmisc)
  
  pixel_info$cluster_label <- umap.res.ref$cluster_label[match(pixel_info$bin, umap.res.ref$bin)]
  
  fc_bin_clean <- HexBin_FC[which(rownames(HexBin_FC) %in% mz_names),]
  
  fc_bin_clean_long <- fc_bin_clean %>%
    rownames_to_column(var = "mz") %>%     
    pivot_longer(
      cols = -mz,                         
      names_to = "HexBin", 
      values_to = "intensity"
    )
  
  fc_bin_clean_long$Label <- pixel_info$cluster_label[match(
    fc_bin_clean_long$HexBin, pixel_info$bin)]
  
  df_pair <- fc_bin_clean_long %>%
    filter(mz %in% c(mz1, mz2)) %>%
    pivot_wider(names_from = mz, values_from = intensity)
  
  df_pair$Label <- factor(df_pair$Label, levels = c("OC", "IC", "Core"))
  
  p1 <- ggplot(df_pair %>%
           filter( !!sym(as.character(mz1)) > 0,
                   !!sym(as.character(mz2)) > 0,
                   is.finite(!!sym(as.character(mz1))),
                   is.finite(!!sym(as.character(mz2))),
                   Label %in% c("OC","IC","Core")),
         aes(x = log2(!!sym(as.character(mz1))), 
             y = log2(!!sym(as.character(mz2))))) +
    geom_density_2d_filled(alpha= 0.95)+
    # scale_fill_brewer(palette= "Spectral", direction = -1) +
    facet_wrap(vars(Label))+
    theme_article() +
    # labs(x = paste("Log2FC of GSH"),
    #      y = paste("Log2FC of CySSG")
    # ) + 
    theme(
      # legend.text = element_text(size = 7),
      # legend.title = element_blank(),
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text  = element_text(size = 8.5)
    )
  
  p2 <- ggplot(df_pair %>%
           filter(!!sym(as.character(mz1)) > 0,
                  !!sym(as.character(mz2)) > 0,
                  is.finite(!!sym(as.character(mz1))),
                  is.finite(!!sym(as.character(mz2))),
                  Label %in% c("OC","IC","Core")),
         aes(x = log2(!!sym(as.character(mz1))), 
             y = log2(!!sym(as.character(mz2))),
             colour = Label)) +
    geom_point(alpha = 0.3, size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    facet_grid(cols = vars(Label)) +
    scale_colour_manual(values = c("OC" = "#692f7c", 
                                   "IC" = "#5C99C5", 
                                   "Core" = "#3C5B8A")) +
    stat_poly_eq(
      aes(label = paste(..rr.label..)),
      formula = y ~ x,
      parse = TRUE,
      label.x.npc = 0.1,   
      label.y.npc = 0.9,   
      size = 3,
      color = "black"
    ) +
    theme_article() +
    # labs(
    #   x = "Log2FC of GSH",
    #   y = "Log2FC of CySSG"
    # ) +
    theme(
      
      legend.title = element_blank(),
      legend.position = "top",
      
    )
  # ggplot(df_pair %>%
  #          filter( !!sym(as.character(mz1)) > 0,
  #                  !!sym(as.character(mz2)) > 0,
  #                  is.finite(!!sym(as.character(mz1))),
  #                  is.finite(!!sym(as.character(mz2))),
  #                  Label %in% c("OC","IC","Core")),
  #        aes(x = log2(!!sym(as.character(mz1))), 
  #            y = log2(!!sym(as.character(mz2))), colour = Label)) +
  #   geom_point(alpha = 0.3, size = 1.5) +
  #   theme_article() +
  #   labs(x = paste("Log2FC of GSH"),
  #        y = paste("Log2FC of CySSG")
  #   ) + 
  #   theme(legend.text = element_text(size = 7),
  #         legend.title = element_blank(),
  #         legend.position = "top",
  #         axis.title = element_text(size = 7),
  #         axis.text  = element_text(size = 5.5) )+
  #   scale_colour_manual(values = c("OC" = "#692f7c", 
  #                                  "IC" = "#5C99C5", 
  #                                  "Core" = "#3C5B8A"))+
  #   geom_smooth(method = "lm", se = FALSE, linewidth = 1) 
  print(p1)
}