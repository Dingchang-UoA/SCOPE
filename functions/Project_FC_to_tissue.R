Project_change_to_tissue <- function(HexBin_FC = HexBin_FC,
                                 mz = 306.0773,
                                 umap.res.ref = umap.res.ref,
                                 pvals_cutoff = F, 
                                 final_pvals = final_pvals,
                                 geom = lens_m, 
                                 type = c("FC", "diff")[1]){
  suppressMessages(library(scales))
  
  if (type == "FC") {
    fc_GSH <- as.data.frame(t(HexBin_FC[which(rownames(HexBin_FC) == mz),]))
    colnames(fc_GSH) <- c("FC")
    fc_GSH$bin <- as.numeric(rownames(fc_GSH))
    
    fc_GSH$bin <- rownames(fc_GSH)
    fc_GSH <- merge(fc_GSH, umap.res.ref %>% 
                      filter(run == "hbnc_a_trim") %>% 
                      select(bin, r.final, z.final, x, y, run),
                    by = c("bin"), all.x = T)
    fc_GSH <- fc_GSH %>% filter(!is.na(run))
    
    outline<- standardize_coordinates(geom)
    outline<- as.data.frame(outline$X)
    colnames(outline) <- c("y", "x")
    outline <- outline %>%
      relocate(x, .before = y)
    outline$group <- 1
    
    fc_vals <- log2(fc_GSH$FC)
    max_abs <- max(abs(fc_vals), na.rm = TRUE)
    limits = c(-3, 3)
    
    
    if (pvals_cutoff == F){
      ggplot(fc_GSH) +
        geom_voronoi(aes(z.final, r.final, fill = fc_vals), outline = outline) +
        coord_equal() +
        theme_article() +
        scale_fill_gradient2(
          low = "#2165aa",        
          mid = "grey90",         
          high = "#a02015",      
          midpoint = 0,
          limits = limits,
          oob = squish             
        ) +
        labs(fill = expression(Log[2]*FC)) +
        theme(
          axis.title = element_blank(),
          axis.text  = element_blank(),
          axis.ticks = element_blank()
        )
      
    } else {
      pvals_flt <- as.data.frame(t(final_pvals[which(rownames(final_pvals) == mz),]))
      pvals_flt$bin <- rownames(pvals_flt)
      colnames(pvals_flt) <- c("pvals", "bin")
      # pvals_flt$FDR <- p.adjust(pvals_flt$pvals)
      
      fc_GSH <- merge(fc_GSH, pvals_flt, by= c("bin"), all.x = T)
      ggplot(fc_GSH) +
        geom_voronoi(
          aes(
            x = z.final,
            y = r.final,
            fill = ifelse(pvals > pvals_cutoff, NA, fc_vals)
          ),
          outline = outline
        ) +
        coord_equal() +
        theme_article() +
        scale_fill_gradient2(
          low = "#2165aa",        
          mid = "grey90",         
          high = "#a02015",      
          midpoint = 0,
          limits = limits,
          oob = squish,
          na.value = "grey90"      
        ) +
        labs(fill = expression(Log[2]*FC)) +
        theme(
          axis.title = element_blank(),
          axis.text  = element_blank(),
          axis.ticks = element_blank()
        )
    }
    
    
    
  } else {
    fc_GSH <- as.data.frame(t(HexBin_FC[which(rownames(HexBin_FC) == mz),]))
    colnames(fc_GSH) <- c("FC")
    fc_GSH$bin <- as.numeric(rownames(fc_GSH))
    
    fc_GSH$bin <- rownames(fc_GSH)
    fc_GSH <- merge(fc_GSH, umap.res.ref %>% 
                      filter(run == "hbnc_a_trim") %>% 
                      select(bin, r.final, z.final, x, y, run),
                    by = c("bin"), all.x = T)
    fc_GSH <- fc_GSH %>% filter(!is.na(run))
    
    outline<- standardize_coordinates(geom)
    outline<- as.data.frame(outline$X)
    colnames(outline) <- c("y", "x")
    outline <- outline %>%
      relocate(x, .before = y)
    outline$group <- 1
    
    ggplot(fc_GSH) +
      geom_voronoi(aes(z.final, r.final, fill = FC), outline = outline) +
      coord_equal() +
      theme_article() +
      labs(fill = "Absolute difference") +
      theme(
        legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank()
      )
    }
  
  
}