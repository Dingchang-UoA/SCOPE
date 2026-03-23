# HiTMaP Compound Scoring Module - Fixed Version

This package contains a fixed version of the HiTMaP `Cpd_spectrum_match_rescore` function and its dependencies.

## Bug Fixes

The original code had the following issues:

1. **Ambiguous column name "formula"**: The error occurred because there were duplicate "formula" columns in different data frames being merged. This has been fixed by using explicit column renaming and specifying the correct columns to join on.

2. **Row count mismatch error (0 vs 1)**: The error occurred when attempting to merge data frames with no matching rows. This has been fixed by adding checks before merges and handling empty data frames gracefully.

3. **Other improvements**:
   - Added more informative error messages
   - Improved handling of column name conflicts
   - Added explicit error handling for potential failure points
   - Fixed potential issues with isotope tables
   - Made column naming more consistent

## Files

The code has been split into modular files to make it easier to maintain:

1. `hitmap-part1-core-functions.R` - Core utilities and basic functions
2. `hitmap-part2-scoring-functions.R` - Functions for spectrum matching and scoring
3. `hitmap-part3-main-function.R` - The main `Cpd_spectrum_match_rescore` function (part 1)
4. `hitmap-part3-main-function-continued.R` - The main function continuation and FDR calculation
5. `hitmap-part4-plotting-functions.R` - Functions for visualization
6. `hitmap-wrapper.R` - Wrapper script to load all components

## Usage

1. Place all files in your working directory
2. Source the wrapper script:

```r
source("hitmap-wrapper.R")
```

3. Call the main function:

```r
result <- Cpd_spectrum_match_rescore(
  cpd_list = your_cpd_list,      # Compound list with formulas
  peaklist = your_peaklist,      # MS peak list data
  wd = "output_directory",       # Where to save results
  SIL_cpd = TRUE,                # Process stable isotope labeled compounds
  adducts = c("M-H", "M+Cl"),    # Adducts to consider
  ppm = 10,                      # PPM tolerance
  score_method = "SQRTP",        # Scoring method
  FDR_cutoff = 0.05              # False discovery rate cutoff
)
```

## Input Data Format

The function expects two main inputs:

1. `cpd_list` - A data frame with the following columns:
   - `formula_mono` or `formula` - Chemical formula
   - `Matched.Form` or `adduct` - Adduct type
   - `Isotype` (optional) - "Normal" or "SIL"
   
2. `peaklist` - A data frame with:
   - `m.z` - m/z values
   - `intensities` - Intensity values

## Debugging Tips

If you encounter issues, you can:

1. Check the message output for error clues
2. Look at intermediate data frames via global assignments:
   - `unique_formula` - The unique formulas being processed
   - `cpd_list_Score` - The scores calculated for each formula
   - `cpd_list_rank` - The ranked compounds

## Dependencies

The code requires several R packages including:

- enviPat
- ggplot2
- stringr
- Rcpp
- dplyr
- Rdisop
- Biostrings
- OrgMassSpecR
- BiocParallel
- data.table
- rcdk
- rcdklibs
- rJava

Use the `check_dependencies()` function to verify all dependencies are installed.
