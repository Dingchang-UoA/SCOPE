remove_duplicates <- function(df, r.e, z.e, r.m, z.m) {
  
  df$distance <- sqrt(df[[r.e]]^2 + df[[z.e]]^2)
  
  df <- df[order(df[[r.m]], df[[z.m]], -df$distance), ]
  
  duplicate_mask <- duplicated(df[, c(r.m, z.m)])
  
  if (any(duplicate_mask)) {
    message("Duplicates found and removed.")
  } else {
    message("No duplicates found.")
  }
  
  # Remove duplicates, keeping the row with max distance
  df <- df[!duplicate_mask, ]
  
  # Drop the temporary distance column
  df$distance <- NULL
  
  return(df)
}
