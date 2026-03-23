Plot_Jitter <- function(data,
                        label_mz_names,
                        mz_column = "m.z",
                        label_column = "metabolite_name",
                        log2FC_column = "log2FC",
                        region_column = "Region",
                        change_column = "change") {
  
  suppressMessages(suppressWarnings(library(ggplot2)))
  suppressMessages(suppressWarnings(library(ggrepel)))
  
  # Add a column for metabolite names to label
  data[[label_column]] <- NA
  
  for (metabolite in names(label_mz_names)) {
    mz_val <- label_mz_names[[metabolite]]
    match_idx <- which(round(as.numeric(data[[mz_column]]), 4) == round(mz_val, 4))
    
    # Only label points that are not "Not Sig"
    match_idx <- match_idx[data[[change_column]][match_idx] != "Not Sig"]
    
    if (length(match_idx) > 0) {
      data[[label_column]][match_idx] <- metabolite
    }
  }
  
  # Create plot
  p <- ggplot(data, aes_string(x = region_column, y = log2FC_column)) +
    geom_jitter(aes_string(color = change_column,
                           size = paste0("abs(", log2FC_column, ")")),
                alpha = 0.45,
                position = position_jitter(width = 0.3)) +
    
    # Add labels with ggrepel (connected lines, no gap)
    geom_text_repel(
      data = subset(data, !is.na(data[[label_column]]) & data[[change_column]] != "Not Sig"),
      aes_string(label = label_column),
      size = 2.54,
      na.rm = TRUE,
      max.overlaps = 10,
      box.padding = 0.4,
      min.segment.length = 0
    ) +
    
    # Customize color and size
    scale_color_manual(values = c("Decreased" = "#8a3684",
                                  "Increased" = "#F7be07",
                                  "Not Sig" = "gray")) +
    scale_size(range = c(2, 5)) +
    
    # Add thresholds
    geom_hline(yintercept = c(-0.5, 0.5),
               linetype = "dashed",
               color = "#999999") +
    
    # Style
    theme_article(base_size = 14) +
    theme(
      legend.position = "top",
      legend.title = element_blank()
    ) +
    guides(
      size = "none"  # Remove size legend
    )
  
  return(p)
}