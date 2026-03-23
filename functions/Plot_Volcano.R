Plot_Volcano <- function(data = results, 
                         pval_cutoff = 0.05, 
                         FC_cutoff = 0.5, 
                         label_mz_names = NULL,
                         label_only_sig = TRUE) {
  
  results <- data
  results$log10_p <- -log10(results$p.value)
  results$log2FC <- log2(results$FC)
  results <- results[is.finite(results$log2FC), ]
  
  results$Regulation <- with(results, 
                             ifelse(log2FC < -1 * FC_cutoff & p.value < pval_cutoff, "Decreased",
                                    ifelse(log2FC > FC_cutoff & p.value < pval_cutoff,
                                           "Increased", "Not Significant")))
  
  # Set factor levels for consistent ordering in legend
  results$Regulation <- factor(results$Regulation, levels = c("Decreased", "Not Significant", "Increased"))
  
  # Add metabolite labels if requested
  results$metabolite_name <- NA
  if (!is.null(label_mz_names)) {
    for (metabolite in names(label_mz_names)) {
      mz_val <- label_mz_names[[metabolite]]
      match_idx <- which(
        round(as.numeric(results$m.z), 4) == round(mz_val, 4) &
          (!label_only_sig | (results$p.value < pval_cutoff & abs(results$log2FC) > FC_cutoff))
      )
      if (length(match_idx) > 0) {
        results$metabolite_name[match_idx] <- metabolite
      }
    }
  }
  
  counts <- table(results$Regulation)
  custom_labels <- c(
    paste("Decreased (", counts["Decreased"], ")", sep = ""),
    paste("Not Significant (", counts["Not Significant"], ")", sep = ""),
    paste("Increased (", counts["Increased"], ")", sep = "")
  )
  
  p <- ggplot(results, aes(x = log2FC, y = log10_p)) +
    geom_point(aes(color = Regulation), alpha = 0.7, size = 3) +
    
    # Decreased labels (to the left)
    geom_text_repel(
      data = subset(results, !is.na(metabolite_name) & Regulation == "Decreased"),
      aes(label = metabolite_name),
      size = 4,
      nudge_x = -0.5,
      na.rm = TRUE,
      max.overlaps = 5,
      box.padding = 0.5,
      min.segment.length = 0
    ) +
    
    # Increased labels (to the right)
    geom_text_repel(
      data = subset(results, !is.na(metabolite_name) & Regulation == "Increased"),
      aes(label = metabolite_name),
      size = 4,
      nudge_x = 0.5,
      na.rm = TRUE,
      max.overlaps = 5,
      box.padding = 0.5,
      min.segment.length = 0
    ) +
    
    scale_color_manual(values = c("Decreased" = "#8a3684",
                                  "Increased" = "#F7be07", 
                                  "Not Significant" = "gray"),
                       labels = custom_labels) +
    theme_minimal() +
    labs(x = "Log2(FC)", y = "-Log10(p-value)") +
    geom_vline(xintercept = c(-1 * FC_cutoff, FC_cutoff), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black") +
    theme_article() +
    theme(legend.position = "top", legend.title = element_blank())
  
  return(p)
}
