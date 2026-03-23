HexBin_FC_for_plot <- function (HexBin_FC = HexBin_FC, mz = mz_names, mz_tol = 30, umap_res, reference = "hbnc"){
  label_mz_df <- enframe(mz, name = "compound", value = "mz")
  
  faceted_fc_data <- pmap_dfr(label_mz_df, function(compound, mz) {
    cpd_idx <- which(
      as.numeric(rownames(HexBin_FC)) %in% match_mz(mz, as.numeric(rownames(HexBin_FC)), tol = mz_tol, unit = "ppm")
    )
    
    if (length(cpd_idx) == 0) return(NULL)  
    
    cpd_fc <- data.frame(
      bin = colnames(HexBin_FC),
      fc = unlist(HexBin_FC[cpd_idx, ]),
      compound = compound
    )
    
    umap_res %>%
      filter(treatment == reference) %>%
      inner_join(cpd_fc, by = "bin")
  })
  
  faceted_fc_data$fc <- as.numeric(faceted_fc_data$fc )
  
  faceted_fc_data_clean <- faceted_fc_data %>%
    filter(!is.na(fc), is.finite(log2(fc)))
  
  return(faceted_fc_data_clean)
}