HexBin_cor <- function (hexbin_df_all,
                        mz_value = 306.0773,
                        run1,
                        run2,
                        runs = c("hbnc_a_trim",
                          "hbnc_b_trim",
                          "hbnc_c_trim")){
  
  
  col_info_new <- tibble(
    colname = colnames(hexbin_df_all),
    run  = str_extract(colnames(hexbin_df_all), ".*(?=@)"),
    hexID   = str_extract(colnames(hexbin_df_all), "(?<=@)\\d+") %>% as.numeric()
  )
  
  
  
  df_long_new <- hexbin_df_all %>%
    tibble::rownames_to_column("mz") %>%
    pivot_longer( -mz,
                  names_to = "colname",
                  values_to = "intensity"
    )  %>%
    left_join(col_info_new, by = "colname") %>%
    filter(mz == as.character(mz_value),
           run %in% runs)

  df_pair_new <- df_long_new %>% 
    filter(run %in% c(run1, run2))
  
  df_mz_new <- df_pair_new %>%
    filter(mz == mz_value) %>%
    select(run, hexID, intensity)
  
  df_wide_new <- df_mz_new %>%
    pivot_wider(names_from = run, 
                values_from = intensity)
  
  df_wide_new[is.na(df_wide_new)] <- 0
  cor(df_wide_new[[run1]], df_wide_new[[run2]])
}