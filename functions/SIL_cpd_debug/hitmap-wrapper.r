### HiTMaP Compound Scoring Module - Wrapper Script

# Load all necessary parts
source("hitmap-part1-core-functions.R")
source("hitmap-part2-scoring-functions.R")
source("hitmap-part3-main-function.R")
source("hitmap-part4-plotting-functions.R")

# Print usage information
cat("HiTMaP Compound Scoring Module loaded successfully.\n\n")
cat("The main functionality is provided by the Cpd_spectrum_match_rescore() function.\n")
cat("Usage example:\n\n")
cat("result <- Cpd_spectrum_match_rescore(\n")
cat("  cpd_list = your_cpd_list,         # Compound list with formulas\n")
cat("  peaklist = your_peaklist,         # MS peak list data\n")
cat("  wd = 'output_directory',          # Where to save results\n")
cat("  SIL_cpd = TRUE,                   # Process stable isotope labeled compounds\n")
cat("  adducts = c('M-H', 'M+Cl'),       # Adducts to consider\n")
cat("  ppm = 10,                         # PPM tolerance\n")
cat("  score_method = 'SQRTP',           # Scoring method\n")
cat("  FDR_cutoff = 0.05                 # False discovery rate cutoff\n")
cat(")\n\n")
cat("The function has been fixed to handle ambiguous column names and to properly join data frames.\n")

# Function to check if all dependencies are installed
check_dependencies <- function() {
  required_packages <- c(
    "enviPat", "ggplot2", "stringr", "Rcpp", "dplyr", "Rdisop", 
    "Biostrings", "OrgMassSpecR", "BiocParallel", "data.table", 
    "rcdk", "rcdklibs", "rJava", "grid", "matter", "OneR", "zoo"
  )
  
  missing <- character(0)
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }
  
  if (length(missing) > 0) {
    cat("The following packages are required but not installed:\n")
    cat(paste(" -", missing), sep = "\n")
    cat("\nPlease install these packages for full functionality.\n")
    cat("You can install them with:\n")
    cat("BiocManager::install(c(", paste(shQuote(missing), collapse = ", "), "))\n")
    return(FALSE)
  } else {
    cat("All dependencies are installed.\n")
    return(TRUE)
  }
}

# Check dependencies when this script is sourced (but don't fail if missing)
tryCatch({
  check_dependencies()
}, error = function(e) {
  cat("Error checking dependencies:", e$message, "\n")
  cat("Some required packages may be missing.\n")
})
