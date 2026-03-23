load_csv <- function(file_path) {
  # Read CSV
  df <- readr::read_csv(file_path)
  
  # Convert to data frame
  df <- as.data.frame(df)
  
  # Set rownames from the first column
  rownames(df) <- df[, 1]
  
  # Drop the first column
  df <- df[, -1, drop = FALSE]
  
  return(df)
}
