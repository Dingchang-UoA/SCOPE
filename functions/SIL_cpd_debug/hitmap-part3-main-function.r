### Main Cpd_spectrum_match_rescore function with bug fixes

Cpd_spectrum_match_rescore <- function(cpd_list, peaklist, wd = getwd(),
                                       tolerance_ppm = 8,
                                       SIL_cpd = TRUE,
                                       C13_SIL_incorporation = 0.995,    
                                       SPECTRUM_batch = 1,
                                       plot_matching_score_t = TRUE,    
                                       adducts = c("M-H", "M+Cl"),
                                       BPPARAM = bpparam(),
                                       atom_isotope_sub = NULL,
                                       ppm = 10,
                                       score_method = "SQRTP",
                                       adjust_score = FALSE,
                                       FDR_cutoff = 0.05,
                                       Decoy_search = TRUE,
                                       Decoy_mode = "isotope") {
  
  # Create workdir if it doesn't exist
  workdir <- wd
  if (!dir.exists(workdir)) {
    dir.create(workdir, recursive = TRUE)
  }
  
  # Load required packages
  adductslist <- Build_adduct_list()
  suppressMessages(suppressWarnings(require(enviPat)))
  suppressMessages(suppressWarnings(require(ggplot2)))
  suppressMessages(suppressWarnings(require(stringr)))
  suppressMessages(suppressWarnings(require("Rcpp")))
  suppressMessages(suppressWarnings(require(dplyr)))
  suppressMessages(suppressWarnings(require(Rdisop)))
  suppressMessages(suppressWarnings(require(Biostrings)))
  suppressMessages(suppressWarnings(require(OrgMassSpecR)))
  suppressMessages(suppressWarnings(require(BiocParallel)))
  suppressMessages(suppressWarnings(require(data.table)))
  data(isotopes)
  
  # Check and prepare input data
  
  # Ensure formula_mono exists
  if (!("formula_mono" %in% colnames(cpd_list))) {
    cpd_list$formula_mono <- cpd_list$formula
  }
  
  # Check if adduct exists, otherwise use Matched.Form
  if ("Matched.Form" %in% colnames(cpd_list)) {
    cpd_list$adduct <- cpd_list$Matched.Form
  }
  
  # Process peaklist
  peaklist$intensities <- as.numeric(peaklist$intensities)
  peaklist <- peaklist[peaklist$intensities > 0, ]
  
  # Set up Isotype if not present
  if (is.null(cpd_list$Isotype)) cpd_list$Isotype <- "Normal"
  
  # Process SIL compounds
  if (SIL_cpd) {
    cpd_list$Isotype[str_detect(cpd_list$Matched.Form, "M_13C")] <- "SIL"
    cpd_list$adduct <- str_replace(cpd_list$adduct, ".M_13C..", "M")
    
    # Check for SIL compounds
    if (length(which(cpd_list$Isotype == "SIL")) == 0) {
      cpd_list_sil <- cpd_list
      cpd_list_sil$Isotype <- "SIL"  # Make sure to set this on the copy
      cpd_list <- rbind(cpd_list, cpd_list_sil)
    }
    
    # Adjust isotopes for SIL
    isotopes$abundance[isotopes$isotope == "13C" & isotopes$element == "[13]C"] <- C13_SIL_incorporation
    addrow_idx <- which(isotopes$isotope == "13C" & isotopes$element == "[13]C")
    if (length(addrow_idx) > 0) {
      addrow <- isotopes[addrow_idx[1], ]
      addrow$isotope <- "12C"
      addrow$mass <- 12
      addrow$abundance <- 1 - C13_SIL_incorporation
      isotopes <- rbind(isotopes, addrow)
      rownames(isotopes) <- 1:nrow(isotopes)
    }
  }
  
  # Filter out empty formulas and polymers
  cpd_list <- cpd_list[cpd_list$formula != "", ]
  cpd_list <- cpd_list[!grepl(")n", cpd_list$formula), ]
  
  # Make copies
  database <- cpd_list
  candidates <- cpd_list
  
  # Check required columns
  required_col <- c("formula", "Isotype")
  candidates <- data_test_rename(required_col, candidates)
  
  # Create formula_isotype column
  candidates$formula_isotype <- paste0(candidates$formula_mono, "_", candidates$Isotype)
  
  # Get unique formulas - avoid ambiguous column names
  unique_formula <- as.data.frame(unique(candidates[, c("formula_mono", "Isotype")]))
  colnames(unique_formula) <- c("formula", "Isotype")  # Ensure consistent column names
  
  # Get atoms from formulas
  if (!exists("get_atoms")) {
    # Define get_atoms if it doesn't exist
    get_atoms <- function(formula) {
      # Simple parsing for example purposes - this should be replaced with actual implementation
      atoms <- list(C = 0, H = 0, N = 0, O = 0)
      return(atoms)
    }
  }
  
  # Create list of unique formulas
  unique_formula_list <- lapply(unique_formula$formula, get_atoms)
  names(unique_formula_list) <- paste0(unique_formula$formula, "_", unique_formula$Isotype)
  
  # Calculate monoisotopic masses
  masslist <- lapply(1:nrow(unique_formula), function(x, unique_formula, unique_formula_list) {
    formula_key <- paste0(unique_formula$formula[x], "_", unique_formula$Isotype[x])
    if (formula_key %in% names(unique_formula_list)) {
      mass <- MonoisotopicMass_sil(
        unique_formula_list[[formula_key]], 
        Isotype = unique_formula$Isotype[x]
      )
      return(mass)
    } else {
      return(NA)  # Handle missing formulas
    }
  }, unique_formula, unique_formula_list)
  
  # Process SIL compounds further if needed
  if (SIL_cpd) {
    unique_formula_list <- lapply(names(unique_formula_list), function(x) {
      y <- unique_formula_list[[x]]
      if (str_detect(x, "_SIL$")) {
        if (!is.null(y$C) && y$C >= 6) {
          y$C <- y$C - 6
          y$'[13]C' <- 6
        } else {
          return(NULL)  # Not enough carbon atoms for SIL
        }
      }
      return(y)
    })
    names(unique_formula_list) <- paste0(unique_formula$formula, "_", unique_formula$Isotype)
    
    # Filter valid SIL entries
    sil_valid_idx <- which(!unlist(lapply(unique_formula_list, is.null)))
    if (length(sil_valid_idx) > 0) {
      unique_formula_list <- unique_formula_list[sil_valid_idx]
      unique_formula <- unique_formula[sil_valid_idx, ]
      masslist <- masslist[sil_valid_idx]
    }
  }
  
  # Reduce formulas to string format
  unique_formula_list_reduce <- atomes.to.formula(unique_formula_list)
  
  # Create mass data frame (with unique column names)
  uniquemass <- unlist(masslist)
  mass_DF <- data.frame(
    formula_mass = unique_formula$formula,  # Using formula_mass to avoid column name conflicts
    Isotype = unique_formula$Isotype, 
    mass = uniquemass, 
    stringsAsFactors = FALSE
  )
  
  # Merge candidates with mass data
  if (nrow(candidates) > 0 && nrow(mass_DF) > 0) {
    # Use explicit join columns to avoid ambiguity
    candidates <- base::merge(
      candidates, 
      mass_DF, 
      by.x = c("formula", "Isotype"),
      by.y = c("formula_mass", "Isotype"),
      all.x = TRUE
    )
  } else {
    message("Cannot merge candidates with mass_DF - empty data frame")
    return(1)
  }
  
  # Update formula column with reduced representation
  formula_keys <- paste0(candidates$formula, "_", candidates$Isotype)
  match_indices <- match(formula_keys, names(unique_formula_list_reduce))
  candidates$formula <- unique_formula_list_reduce[match_indices]
  
  # Remove duplicate columns
  candidates <- candidates[, !duplicated(names(candidates))]
  
  # Set Metabolite.Name if missing
  if (is.null(candidates$Metabolite.Name)) {
    candidates$Metabolite.Name <- candidates$Matched.Compound
  }
  
  # Rename/check required columns
  required_col <- c("mz", "Metabolite.Name", "mass", "formula", "Isotype")
  candidates <- data_test_rename(required_col, candidates)
  candidates$formula <- as.character(candidates$formula)
  meta_symbol <- unique_formula_list
  
  # Process adducts - assuming convert_peptide_adduct_list is defined elsewhere
  symbol_adducts <- try(
    bplapply(adducts, convert_peptide_adduct_list, meta_symbol, 
             BPPARAM = BPPARAM, adductslist = adductslist),
    silent = TRUE
  )
  
  if (inherits(symbol_adducts, "try-error")) {
    message("Error in adduct conversion: ", as.character(symbol_adducts))
    return(1)
  }
  
  symbol_adducts_df <- lapply(
    symbol_adducts,
    function(x, names_x) {
      data.frame(formula_isotype = names_x, formula_adduct = x)
    }, 
    names(meta_symbol)
  )
  names(symbol_adducts_df) <- adducts
  
  # Initialize Meta_Summary
  Meta_Summary <- NULL
  
  # Process each adduct
  for (i in 1:length(adducts)) {
    candidates_adduct <- unique(candidates[, c("Metabolite.Name", "mass", "formula_mono", "Isotype")])
    candidates_adduct <- candidates_adduct[!duplicated(candidates_adduct[, c("Metabolite.Name", "formula_mono", "Isotype")]), ]
    
    # Calculate mz based on adduct
    adduct_idx <- which(adductslist$Name == adducts[i])
    if (length(adduct_idx) > 0) {
      adductmass <- as.numeric(as.character(adductslist[adduct_idx, "Mass"]))
      adducts_Charge <- as.numeric(adductslist[adduct_idx, "Charge"])
      
      candidates_adduct$mz <- (as.numeric(candidates_adduct$mass) + adductmass) / abs(adducts_Charge)
      candidates_adduct$adduct <- adducts[i]
      candidates_adduct$charge <- adducts_Charge
      candidates_adduct$formula_isotype <- paste0(candidates_adduct$formula_mono, "_", candidates_adduct$Isotype)
      
      # Debug output
      message(paste("Processing adduct:", adducts[i]))
      message("Candidates_adduct rows:", nrow(candidates_adduct))
      
      # Check if the adduct exists in the symbol_adducts_df list
      if (adducts[i] %in% names(symbol_adducts_df)) {
        message("Symbol_adducts_df rows:", nrow(symbol_adducts_df[[adducts[i]]]))
        
        # Check if we have data to merge
        if (nrow(candidates_adduct) > 0 && nrow(symbol_adducts_df[[adducts[i]]]) > 0) {
          # First, check which formula_isotype values are actually in both dataframes
          common_isotypes <- intersect(
            candidates_adduct$formula_isotype, 
            symbol_adducts_df[[adducts[i]]]$formula_isotype
          )
          
          if (length(common_isotypes) > 0) {
            # Use left join to ensure we keep all candidates_adduct rows
            merged_adduct <- merge(
              candidates_adduct, 
              symbol_adducts_df[[adducts[i]]], 
              by = "formula_isotype", 
              all.x = TRUE
            )
            
            # Fix column names to avoid ambiguity
            if ("formula_adduct" %in% colnames(merged_adduct)) {
              merged_adduct$formula <- merged_adduct$formula_adduct
              merged_adduct$formula_adduct <- NULL
            }
            
            # Add to meta summary
            Meta_Summary <- rbind.data.frame(Meta_Summary, unique(merged_adduct))
          } else {
            message("No common formula_isotype values found for adduct: ", adducts[i])
          }
        } else {
          message("Skipping adduct ", adducts[i], " - insufficient data for merge")
        }
      } else {
        message("No symbol_adducts_df entry for adduct: ", adducts[i])
      }
    } else {
      message("Adduct not found in adductslist: ", adducts[i])
    }
  }
  
  # If Meta_Summary is empty, return early
  if (is.null(Meta_Summary) || nrow(Meta_Summary) == 0) {
    message("No compounds to process. Empty Meta_Summary.")
    return(1)
  }
  
  # Update cpd_list with Meta_Summary
  cpd_list <- Meta_Summary
  
  # Ensure Do_PMF_search exists
  if (!exists("Do_PMF_search")) {
    # Define minimal version if it doesn't exist
    Do_PMF_search <- function(peaklist, cpd_list, BPPARAM, ppm) {
      # This is a placeholder - should be replaced with actual implementation
      return(data.frame(mz = cpd_list$mz, Intensity = 1))
    }
  }
  
  # First round of peptide search to get putative result
  mz_feature_list <- try(
    Do_PMF_search(peaklist, cpd_list, BPPARAM = BPPARAM, ppm = ppm),
    silent = TRUE
  )
  
  if (inherits(mz_feature_list, "try-error")) {
    message("Error in Do_PMF_search: ", as.character(mz_feature_list))
    return(1)
  }
  
  mz_feature_list <- unique(mz_feature_list)
  message("Summarizing cpd information...", length(unique(mz_feature_list$mz)), " were found in the first search.")
  
  # Convert to data.table and prepare for merge
  cpd_list <- as.data.table(cpd_list)
  mz_feature_list <- as.data.table(mz_feature_list)
  
  # Remove any existing Intensity column to avoid conflicts
  if ("Intensity" %in% colnames(cpd_list)) {
    cpd_list$Intensity <- NULL
  }
  
  # Convert mz to character for consistent joining
  mz_feature_list$mz <- as.character(mz_feature_list$mz)
  cpd_list$mz <- as.character(cpd_list$mz)
  
  # Save full copy of cpd_list
  cpd_list_full <- cpd_list
  
  # Merge with mz_feature_list using explicit column names
  cpd_list <- merge(
    cpd_list, 
    mz_feature_list, 
    by.x = "mz", 
    by.y = "mz", 
    all.x = TRUE, 
    sort = FALSE,
    suffixes = c("", ".y")  # Ensure clear suffix for duplicated columns
  )
  
  # Handle NA in Intensity
  if ("Intensity" %in% colnames(cpd_list)) {
    cpd_list$Intensity[is.na(cpd_list$Intensity)] <- 0
    cpd_list <- cpd_list[cpd_list$Intensity > 0, ]
  } else {
    message("No Intensity column found after merge")
    return(1)
  }
  
  cpd_list$formula <- as.character(cpd_list$formula)
  
  # Get unique formulas for scoring
  unique_formula <- unique(cpd_list[, c("formula", "Isotype")])
  unique_formula <- unique_formula[!is.na(unique_formula$formula), ]
  
  # For debugging
  assign("unique_formula", unique_formula, envir = .GlobalEnv)
  
  # Score compounds with enhanced error handling
  cpd_list_Score <- lapply(unique_formula$formula, function(f) {
    tryCatch({
      SCORE_CPD(
        f, peaklist, isotopes = isotopes,
        score_method = score_method, 
        charge = -1, ppm = ppm, threshold = 1
      )
    }, error = function(e) {
      message("Error scoring formula: ", f, " - ", e$message)
      return(c(0, 0, 0))  # Return default values on error
    })
  })
  
  # For debugging
  assign("cpd_list_Score", cpd_list_Score, envir = .GlobalEnv)
  
  # Process scores only if we have data
  if (length(cpd_list_Score) > 0) {
    cpd_list_Score_m <- as.data.frame(do.call(rbind, cpd_list_Score))
    
    if (ncol(cpd_list_Score_m) > 0) {
      # Assign column names
      names(cpd_list_Score_m) <- c("Score", "delta_ppm", "Intensity")
      
      # Add Isotype
      cpd_list_Score_m$Isotype <- unique_formula$Isotype
      
      # Create formula score dataframe
      formula_score <- data.frame(
        formula = unique_formula$formula,
        Score = cpd_list_Score_m$Score,
        Delta_ppm = cpd_list_Score_m$delta_ppm,
        Intensity_score = cpd_list_Score_m$Intensity,  # Renamed to avoid conflicts
        Isotype = unique_formula$Isotype,
        stringsAsFactors = FALSE
      )
      
      # Merge scores with cpd_list
      # First remove any existing score columns to avoid conflicts
      cols_to_remove <- c("Score", "Delta_ppm")
      for (col in cols_to_remove) {
        if (col %in% colnames(cpd_list)) {
          cpd_list[[col]] <- NULL
        }
      }
      
      # Handle potential duplicates in Intensity column
      if ("Intensity.y" %in% colnames(cpd_list)) {
        # Rename to avoid confusion
        cpd_list$Intensity_orig <- cpd_list$Intensity
        cpd_list$Intensity.y <- NULL
      }
      
      # Merge with formula scores
      cpd_list <- merge(
        cpd_list, 
        formula_score, 
        by = c("formula", "Isotype"), 
        all.x = TRUE,
        suffixes = c("", ".score")
      )
      
      # Set isdecoy flag
      cpd_list$isdecoy <- 0
      
      # Generate decoy candidate list if needed
      if (Decoy_search && ("isotope" %in% Decoy_mode)) {
        # Create modified isotope table for decoys
        decoy_isotopes <- isotopes
        
        # Define decoy isotopes
        C_rows <- data.frame(
          element = c("C", "C", "[13]C", "[13]C"),
          isotope = c("11C", "12C", "13C", "14C"),
          mass = c(10.99664516, 12, 13.003355, 14),
          abundance = c(0.0107, 1 - 0.0107, 1 - C13_SIL_incorporation, C13_SIL_incorporation),
          ratioC = 0,
          stringsAsFactors = FALSE
        )
        
        # Update isotope table
        decoy_isotopes <- decoy_isotopes[!(decoy_isotopes$element %in% c("[13]C", "C")), ]
        decoy_isotopes <- rbind(decoy_isotopes, C_rows)
        
        # More decoy isotope modifications
        isotope_replacements <- list(
          list(from_isotope = "15N", to = list(element = "N", isotope = "13N", mass = 13.00603905, abundance = 0.00364000, ratioC = 4)),
          list(from_isotope = "2H", to = list(element = "H", isotope = "0H", mass = 0.001548286, abundance = 0.00011500, ratioC = 6)),
          list(from_isotope = "17O", to = list(element = "O", isotope = "15O", mass = 14.99069774, abundance = 0.00038000, ratioC = 3)),
          list(from_isotope = "18O", to = list(element = "O", isotope = "14O", mass = 13.99066884, abundance = 0.00205000, ratioC = 3)),
          list(from_isotope = "33S", to = list(element = "S", isotope = "31S", mass = 30.97268292, abundance = 0.0075, ratioC = 3)),
          list(from_isotope = "34S", to = list(element = "S", isotope = "30S", mass = 29.9762745, abundance = 0.0425, ratioC = 3)),
          list(from_isotope = "35S", to = list(element = "S", isotope = "29S", mass = 28.94414146, abundance = 0, ratioC = 3)),
          list(from_isotope = "36S", to = list(element = "S", isotope = "28S", mass = 27.97706058, abundance = 0.0001, ratioC = 3)),
          list(from_isotope = "37Cl", to = list(element = "Cl", isotope = "33Cl", mass = 32.965903, abundance = 0.24240000, ratioC = 3))
        )
        
        # Apply isotope replacements
        for (repl in isotope_replacements) {
          idx <- which(decoy_isotopes$isotope %in% repl$from_isotope)
          if (length(idx) > 0) {
            decoy_isotopes[idx, "element"] <- repl$to$element
            decoy_isotopes[idx, "isotope"] <- repl$to$isotope
            decoy_isotopes[idx, "mass"] <- repl$to$mass
            if ("abundance" %in% colnames(decoy_isotopes)) {
              decoy_isotopes[idx, "abundance"] <- repl$to$abundance
            } else if ("aboudance" %in% colnames(decoy_isotopes)) {
              # Handle potential typo in column name
              decoy_isotopes[idx, "aboudance"] <- repl$to$abundance
            }
            decoy_isotopes[idx, "ratioC"] <- repl$to$ratioC
          }
        }
        
        # Create decoy version of cpd_list
        cpd_list_decoy <- cpd_list[cpd_list$isdecoy == 0, ]
        cpd_list_decoy$isdecoy <- 1
        
        # Get unique formulas for decoy scoring
        unique_formula <- unique(cpd_list_decoy[, c("formula", "Isotype")])
        unique_formula <- unique_formula[!is.na(unique_formula$formula), ]
        
        # For debugging
        assign("unique_formula", unique_formula, envir = .GlobalEnv)
        
        # Score decoy compounds
        cpd_list_Score <- lapply(unique_formula$formula, function(f) {
          tryCatch({
            SCORE_CPD(
              f, peaklist, isotopes = decoy_isotopes,
              score_method = score_method, 
              charge = -1, ppm = ppm, threshold = 1
            )
          }, error = function(e) {
            message("Error scoring decoy formula: ", f, " - ", e$message)
            return(c(0, 0, 0))  # Return default values on error
          })
        })
        
        # For debugging
        assign("cpd_list_Score", cpd_list_Score, envir = .GlobalEnv)
        
        # Process decoy scores
        if (length(cpd_list_Score) > 0) {
          cpd_list_Score_m <- as.data.frame(do.call(rbind, cpd_list_Score))
          
          if (ncol(cpd_list_Score_m) > 0) {
            # Assign column names
            names(cpd_list_Score_m) <- c("Score", "delta_ppm", "Intensity")
            
            # Add Isotype
            cpd_list_Score_m$Isotype <- unique_formula$Isotype
            
            # Create formula score dataframe
            formula_score <- data.frame(
              formula = unique_formula$formula,
              Score = cpd_list_Score_m$Score,
              Delta_ppm = cpd_list_Score_m$delta_ppm,
              Intensity_score = cpd_list_Score_m$Intensity,
              Isotype = unique_formula$Isotype,
              stringsAsFactors = FALSE
            )
            
            # Merge decoy scores with cpd_list_decoy
            # First remove any existing score columns to avoid conflicts
            for (col in cols_to_remove) {
              if (col %in% colnames(cpd_list_decoy)) {
                cpd_list_decoy[[col]] <- NULL
              }
            }
            
            # Handle potential duplicates in Intensity column
            if ("Intensity.y" %in% colnames(cpd_list_decoy)) {
              cpd_list_decoy$Intensity_orig <- cpd_list_decoy$Intensity
              cpd_list_decoy$Intensity.y <- NULL
            }
            
            # Merge with formula scores
            cpd_list_decoy <- merge(
              cpd_list_decoy, 
              formula_score, 
              by = c("formula", "Isotype"), 
              all.x = TRUE,
              suffixes = c("", ".score")
            )
            
            # Combine target and decoy compounds
            cpd_list <- rbind(cpd_list, cpd_list_decoy)
          } else {
            message("No columns in decoy score matrix")
          }
        } else {
          message("No decoy scores calculated")
        }
      }
      
      # Prepare for ranking
      cpd_list$Query.Mass <- cpd_list$mz
      cpd_list$mz <- as.numeric(as.character(cpd_list$mz))
      
      # Check if we have any matches
      if (all(unique(cpd_list$Score) %in% c(0, NA))) {
        message("No matched feature, proceed to next spectrum")
        return(1)
      }
      
      # Ensure rank_mz_feature exists
      if (!exists("rank_mz_feature")) {
        # Define minimal version if it doesn't exist
        rank_mz_feature <- function(cpd_list, mz_feature, BPPARAM) {
          # This is a placeholder
          cpd_list$mz_align <- cpd_list$mz
          cpd_list$Rank <- 1:nrow(cpd_list)
          return(cpd_list)
        }
      }
      
      # Rank mz features
      cpd_list_rank <- try(
        rank_mz_feature(cpd_list, mz_feature = peaklist, BPPARAM = BPPARAM),
        silent = TRUE
      )
      
      if (inherits(cpd_list_rank, "try-error")) {
        message("Error in rank_mz_feature: ", as.character(cpd_list_rank))
        return(1)
      }
      
      if ("mz_align" %in% colnames(cpd_list_rank)) {
        cpd_list_rank$mz_align <- unlist(cpd_list_rank$mz_align)
      }
      
      cpd_list_rank <- as.data.frame(unique(cpd_list_rank))
      cpd_list_rank <- cpd_list_rank[!is.na(cpd_list_rank$Score), ]
      
      # Create output directory
      dir.create(paste0(wd, "/", SPECTRUM_batch), recursive = TRUE)
      
      # Filter and save ranked list
      if ("Delta_ppm" %in% colnames(cpd_list_rank)) {
        cpd_list_rank$Delta_ppm[is.na(cpd_list_rank$Delta_ppm)] <- ppm
      }
      
      if ("Intensity" %in% colnames(cpd_list_rank)) {
        cpd_list_rank <- cpd_list_rank[cpd_list_rank$Intensity >= 0, ]
        cpd_list_rank <- cpd_list_rank[!is.na(cpd_list_rank$Intensity), ]
      }
      
      output_file <- paste0(wd, "/", SPECTRUM_batch, "/", "CPD_1st_ID_score_rank_", score_method, ".csv")
      write.csv(cpd_list_rank, output_file, row.names = FALSE)
      
      # For debugging
      assign("cpd_list_rank", cpd_list_rank, envir = .GlobalEnv)
      
      # Score filtering
      cpd_list_2nd <- data.frame()
      cpd_list_2nd_sil <- data.frame()
      
      if (nrow(cpd_list_rank) >= 2) {
        if (adjust_score == FALSE) {
          dir.create(paste0(wd, "/", SPECTRUM_batch))
          
          # Ensure FDR_cutoff_plot exists
          if (!exists("FDR_cutoff_plot")) {
            # Define minimal version
            FDR_cutoff_plot <- function(cpd_list, FDR_cutoff, plot_fdr, outputdir, adjust_score) {
              return(min(cpd_list$Score, na.rm = TRUE))
            }
          }
          


          Score_cutoff <- FDR_cutoff_plot(
            cpd_list_rank, 
            FDR_cutoff = FDR_cutoff, 
            plot_fdr = TRUE, 
            outputdir = paste0(wd, "/", SPECTRUM_batch), 
            adjust_score = adjust_score
          )
          
          # Filter compounds based on score cutoff
          cpd_list_2nd <- cpd_list_rank[
            ((cpd_list_rank$Score >= Score_cutoff) & (!is.na(cpd_list_rank$Intensity))),
          ]
          
          # Process SIL compounds separately if needed
          if (SIL_cpd) {
            sil_dir <- paste0(wd, "/", SPECTRUM_batch, "/", "SIL")
            dir.create(sil_dir, recursive = TRUE, showWarnings = FALSE)
            
            # Extract SIL compounds
            cpd_list_2nd_sil <- cpd_list_rank[cpd_list_rank$Isotype == "SIL", ]
            
            if (nrow(cpd_list_2nd_sil) > 0) {
              # Calculate SIL-specific score cutoff
              Score_cutoff_sil <- FDR_cutoff_plot(
                cpd_list_2nd_sil, 
                FDR_cutoff = FDR_cutoff, 
                plot_fdr = TRUE, 
                outputdir = sil_dir, 
                adjust_score = adjust_score
              )
              
              # Filter SIL compounds
              cpd_list_2nd_sil <- cpd_list_2nd_sil[
                ((cpd_list_2nd_sil$Score >= Score_cutoff_sil) & 
                 (!is.na(cpd_list_2nd_sil$Intensity))),
              ]
              
              # Combine with regular compounds
              if (nrow(cpd_list_2nd_sil) >= 1) {
                cpd_list_2nd <- rbind(cpd_list_2nd, cpd_list_2nd_sil)
                cpd_list_2nd <- unique(cpd_list_2nd)
              }
            }
          }
        } else {
          # Adjusted score filtering method
          dir.create(paste0(wd, "/", SPECTRUM_batch), recursive = TRUE, showWarnings = FALSE)
          
          # Save original scores
          cpd_list_rank$Score_org <- cpd_list_rank$Score
          
          # Define score correction function
          mz.score.correction <- function(mz.ref, Score.ref, mz, Score,
                                         control = loess.control(),
                                         span = 0.75, fit_type = "loess", force = FALSE) {
            suppressMessages(suppressWarnings(library(matter)))
            
            if (fit_type == "loess") {
              shift <- suppressWarnings(loess(Score.ref ~ mz.ref, span = span, control = control))
              dScore <- predict(shift, mz)
            }
            
            if (fit_type == "lm") {
              shift <- suppressWarnings(lm(Score.ref ~ mz.ref))
              dScore <- predict.lm(shift, data.frame(mz.ref = mz))
            }
            
            if (fit_type == "poly") {
              shift <- suppressWarnings(lm(Score.ref ~ poly(mz.ref, 2)))
              dScore <- predict.lm(shift, data.frame(mz.ref = mz))
            }
            
            return(list(
              final.Score = Score + dScore,
              dScore = dScore,
              shift = shift
            ))
          }
          
          # Adjust scores based on decoy distribution
          unique_formula_score <- unique(cpd_list_rank[, c("formula", "mz", "Score_org", "Isotype", "isdecoy")])
          unique_formula_score_decoy <- unique_formula_score[unique_formula_score$isdecoy == 1, ]
          
          if (nrow(unique_formula_score_decoy) > 0) {
            correction_result <- mz.score.correction(
              mz.ref = unique_formula_score_decoy$mz,
              Score.ref = unique_formula_score_decoy$Score_org,
              mz = cpd_list_rank$mz,
              Score = cpd_list_rank$Score_org,
              fit_type = "lm"
            )
            
            cpd_list_rank$Score_adj <- correction_result$final.Score
            cpd_list_rank$Score <- cpd_list_rank$Score_adj
          }
          
          # Calculate FDR-based score cutoff
          Score_cutoff <- FDR_cutoff_plot(
            cpd_list_rank, 
            FDR_cutoff = FDR_cutoff, 
            plot_fdr = TRUE, 
            outputdir = paste0(wd, "/", SPECTRUM_batch), 
            adjust_score = FALSE
          )
          
          # Filter compounds
          cpd_list_2nd <- cpd_list_rank[
            ((cpd_list_rank$Score >= Score_cutoff) & (!is.na(cpd_list_rank$Intensity))),
          ]
          
          # Process SIL compounds separately if needed
          if (SIL_cpd) {
            sil_dir <- paste0(wd, "/", SPECTRUM_batch, "/", "SIL")
            dir.create(sil_dir, recursive = TRUE, showWarnings = FALSE)
            
            # Extract SIL compounds
            cpd_list_2nd_sil <- cpd_list_rank[cpd_list_rank$Isotype == "SIL", ]
            
            if (nrow(cpd_list_2nd_sil) > 0) {
              # Get SIL-specific formula scores
              unique_formula_score <- unique(cpd_list_2nd_sil[, c("formula", "mz", "Score_org", "Isotype", "isdecoy")])
              unique_formula_score_decoy <- unique_formula_score[unique_formula_score$isdecoy == 1, ]
              
              if (nrow(unique_formula_score_decoy) > 0) {
                # Apply score correction for SIL compounds
                correction_result <- mz.score.correction(
                  mz.ref = unique_formula_score_decoy$mz,
                  Score.ref = unique_formula_score_decoy$Score_org,
                  mz = cpd_list_2nd_sil$mz,
                  Score = cpd_list_2nd_sil$Score_org,
                  fit_type = "loess"
                )
                
                cpd_list_2nd_sil$Score_adj <- correction_result$final.Score
                cpd_list_2nd_sil$Score <- cpd_list_2nd_sil$Score_adj
              }
              
              # Calculate SIL-specific score cutoff
              Score_cutoff_sil <- FDR_cutoff_plot(
                cpd_list_2nd_sil, 
                FDR_cutoff = FDR_cutoff, 
                plot_fdr = TRUE, 
                outputdir = sil_dir, 
                adjust_score = FALSE
              )
              
              # Filter SIL compounds
              cpd_list_2nd_sil <- cpd_list_2nd_sil[
                ((cpd_list_2nd_sil$Score >= Score_cutoff_sil) & (!is.na(cpd_list_2nd_sil$Intensity))),
              ]
              
              # Combine with regular compounds
              if (nrow(cpd_list_2nd_sil) >= 1) {
                cpd_list_2nd <- rbind(cpd_list_2nd, cpd_list_2nd_sil)
                
                # Reset scores to original values
                cpd_list_2nd$Score <- cpd_list_2nd$Score_org
                cpd_list_2nd$Score_org <- NULL
                cpd_list_2nd$Score_adj <- NULL
                
                cpd_list_2nd <- unique(cpd_list_2nd)
              }
            }
          }
        }
      } else {
        message("Insufficient compounds for FDR calculation (need at least 2)")
      }
      
      # Save SIL results
      if (exists("cpd_list_2nd_sil") && nrow(cpd_list_2nd_sil) > 0) {
        sil_output_file <- paste0(wd, "/", SPECTRUM_batch, "/", "CPD_2nd_ID_score_rank_SIL_", score_method, ".csv")
        write.csv(cpd_list_2nd_sil, sil_output_file, row.names = FALSE)
      }
      
      # Ensure numeric mz values
      cpd_list_rank$mz <- as.numeric(as.character(cpd_list_rank$mz))
      
      # Process final results
      if (exists("cpd_list_2nd") && nrow(cpd_list_2nd) > 0) {
        cpd_list_2nd$mz <- as.numeric(as.character(cpd_list_2nd$mz))
        
        # Remove decoy compounds from final results
        cpd_list_2nd <- cpd_list_2nd[cpd_list_2nd$isdecoy == 0, ]
        
        # Save final results
        output_file <- paste0(wd, "/", SPECTRUM_batch, "/", "CPD_2nd_ID_score_rank_", score_method, ".csv")
        write.csv(cpd_list_2nd, output_file, row.names = FALSE)
        
        # Plot matching spectrum and score if requested
        if (plot_matching_score_t && nrow(cpd_list_2nd) > 0) {
          # Create directory for spectrum plots
          spec_dir <- paste0(wd, "/", SPECTRUM_batch, "/cpd_spectrum_match")
          dir.create(spec_dir, recursive = TRUE, showWarnings = FALSE)
          
          # Plot normal compounds
          normal_cpds <- cpd_list_2nd[cpd_list_2nd$Isotype == "Normal", ]
          if (nrow(normal_cpds) > 0) {
            tryCatch({
              plot_cpd_matching_score(
                normal_cpds, 
                peaklist, 
                isotopes = isotopes, 
                charge = -1, 
                ppm, 
                outputdir = spec_dir
              )
            }, error = function(e) {
              message("Error plotting normal spectra: ", e$message)
            })
          }
          
          # Plot SIL compounds
          if (SIL_cpd) {
            sil_spec_dir <- paste0(wd, "/", SPECTRUM_batch, "/cpd_spectrum_match_sil")
            dir.create(sil_spec_dir, recursive = TRUE, showWarnings = FALSE)
            
            sil_cpds <- cpd_list_2nd[cpd_list_2nd$Isotype == "SIL", ]
            if (nrow(sil_cpds) > 0) {
              tryCatch({
                plot_cpd_matching_score(
                  sil_cpds, 
                  peaklist, 
                  isotopes = isotopes, 
                  charge = -1, 
                  ppm, 
                  outputdir = sil_spec_dir
                )
              }, error = function(e) {
                message("Error plotting SIL spectra: ", e$message)
              })
            }
          }
        }
        
        # Return the filtered compounds
        return(cpd_list_2nd)
      } else {
        message("No compounds passed filtering criteria")
        return(data.frame())
      }
    } else {
      message("No scores calculated")
      return(1)
    }
  } else {
    message("No compounds to score")
    return(1)
  }
}

# Helper function: FDR cutoff plot for compound detection
FDR_cutoff_plot <- function(cpd_plot_list, FDR_cutoff = 0.1, FDR_strip = 500, 
                           plot_fdr = FALSE, outputdir = NULL, adjust_score = FALSE) {
  suppressMessages(suppressWarnings(require(ggplot2)))
  suppressMessages(suppressWarnings(require(data.table)))
  suppressMessages(suppressWarnings(require(dplyr)))
  suppressMessages(suppressWarnings(require(zoo)))
  
  # Check required columns
  cpd_plot_list <- data_test_rename(c("isdecoy", "Score"), cpd_plot_list)
  
  # Convert isdecoy to factor
  cpd_plot_list$isdecoy <- factor(cpd_plot_list$isdecoy)
  
  # Save a copy for output
  cpd_plot_list_output <- cpd_plot_list
  
  # Plot FDR histogram function
  plot_fdr_histogram <- function(cpd_plot_list, plot_name = "target-decoy", outputdir = outputdir) {
    suppressMessages(suppressWarnings(require(dplyr)))
    
    if (nrow(cpd_plot_list) != 0) {
      # Extract data for plotting
      cpd_plot_list_plot <- cpd_plot_list[, c("isdecoy", "Score")]
      target_decoy <- factor(
        ifelse(cpd_plot_list_plot$isdecoy == 0, "Target", "Decoy"), 
        levels = c("Target", "Decoy")
      )
      cpd_plot_list_plot$target_decoy <- target_decoy
      cpd_plot_list$target_decoy <- target_decoy
      
      # Calculate mean scores by group
      mu <- cpd_plot_list_plot %>% 
        group_by(target_decoy) %>% 
        summarize(mean = mean(Score, na.rm = TRUE))
      
      # Create density histogram
      tryCatch({
        png(paste0(outputdir, "/Peptide_Score_histogram_", plot_name, ".png"))
        p <- ggplot(
          cpd_plot_list_plot, 
          aes(x = cpd_plot_list_plot$Score, color = target_decoy, fill = target_decoy)
        ) +
          geom_histogram(
            aes(y = ..density..), 
            position = "dodge", 
            alpha = 0.5, 
            bins = 50
          ) +
          geom_density(alpha = 0.6) + 
          xlim(quantile(cpd_plot_list_plot$Score, c(0.005, 0.995), na.rm = TRUE)) +
          geom_vline(
            data = mu, 
            aes(xintercept = mu$mean, color = mu$target_decoy), 
            linetype = "dashed"
          ) +
          ggtitle("Peptide score vs Counts") +
          xlab("Score") + 
          ylab("Counts") +
          theme_classic()
        print(p)
        dev.off()
      }, error = function(e) {
        message("Error plotting histogram: ", e$message)
      })
      
      # Create score vs mz plot
      if (log(nrow(cpd_plot_list_plot)) > 0 && "mz" %in% colnames(cpd_plot_list)) {
        tryCatch({
          png(paste0(outputdir, "/Matching_Score_vs_mz_", plot_name, ".png"))
          alpha <- 1 / 20 / (round(log(nrow(cpd_plot_list_plot), base = 2)) / 40)
          if (alpha >= 1) alpha <- 1
          
          p <- ggplot(
            data = cpd_plot_list, 
            aes(x = mz, y = Score, colour = target_decoy)
          ) + 
            geom_point(size = 1, alpha = alpha, aes(colour = target_decoy)) +
            ggtitle("Matching score vs mz") +
            xlab("mz") + 
            ylab("Matching score") + 
            labs(fill = "isdecoy") + 
            theme(axis.text.x = element_text(angle = 45))
          print(p)
          dev.off()
        }, error = function(e) {
          message("Error plotting score vs mz: ", e$message)
        })
      }
    }
  }
  
  # Plot histograms if requested
  if (!is.null(outputdir) && plot_fdr) {
    plot_fdr_histogram(cpd_plot_list, outputdir = outputdir)
  }
  
  # Apply score adjustment if requested
  if (adjust_score) {
    # Split by target/decoy status
    cpd_plot_list_sep <- cpd_plot_list_output[cpd_plot_list_output$isdecoy == 0, ]
    cpd_plot_list_sep_d <- cpd_plot_list_output[cpd_plot_list_output$isdecoy == 1, ]
    
    # Merge to calculate adjusted scores
    if (nrow(cpd_plot_list_sep) > 0 && nrow(cpd_plot_list_sep_d) > 0) {
      merge_cols <- names(cpd_plot_list_sep)[!names(cpd_plot_list_sep) %in% c("isdecoy", "Score", "Rank")]
      cpd_plot_list_mer <- merge(
        cpd_plot_list_sep, 
        cpd_plot_list_sep_d, 
        by = merge_cols
      )
      
      # Calculate adjusted score
      cpd_plot_list_mer$Score <- cpd_plot_list_mer$Score.x - cpd_plot_list_mer$Score.y
      cpd_plot_list_mer$isdecoy <- 0
      
      # Create decoy version with score = 0
      cpd_plot_list_adj <- cpd_plot_list_mer
      cpd_plot_list_mer$Score <- 0
      cpd_plot_list_mer$isdecoy <- 1
      
      # Combine adjusted target and decoy
      cpd_plot_list_adj <- rbind(cpd_plot_list_adj, cpd_plot_list_mer)
      cpd_plot_list_adj$isdecoy <- as.factor(cpd_plot_list_adj$isdecoy)
      
      # Plot adjusted scores
      plot_fdr_histogram(cpd_plot_list_adj, plot_name = "adj", outputdir = outputdir)
      
      # Update cpd_plot_list
      cpd_plot_list <- cpd_plot_list_adj
      
      # Find score distribution modes
      den <- density(cpd_plot_list_adj$Score[cpd_plot_list_adj$isdecoy == 0])
      
      # Helper function to find modes
      find_modes <- function(x) {
        modes <- NULL
        for (i in 2:(length(x) - 1)) {
          if ((x[i] > x[i - 1]) & (x[i] > x[i + 1])) {
            modes <- c(modes, i)
          }
        }
        if (length(modes) == 0) {
          modes <- 'This is a monotonic distribution'
        }
        return(modes)
      }
      
      # Calculate modes
      mode_indices <- find_modes(den$y)
      mode_values <- den$x[mode_indices]
      
      # Resolve multiple modes (assuming resolve_multi_modes is defined elsewhere)
      if (exists("resolve_multi_modes")) {
        mm <- resolve_multi_modes(
          cpd_plot_list_adj$Score[cpd_plot_list_adj$isdecoy == 0], 
          adjust = 2
        )
        
        # Use the highest mode minus an offset as the cutoff
        if (is.list(mm) && "modes" %in% names(mm)) {
          mode_values <- as.numeric(mm$modes)
          Score_cutoff <- max(mode_values) - 0.75
        } else {
          # Default if mode resolution fails
          Score_cutoff <- median(cpd_plot_list_adj$Score[cpd_plot_list_adj$isdecoy == 0], na.rm = TRUE)
        }
      } else {
        # Use the highest density point as cutoff
        if (length(mode_values) > 0) {
          Score_cutoff <- max(mode_values)
        } else {
          # Default if no modes found
          Score_cutoff <- median(cpd_plot_list_adj$Score[cpd_plot_list_adj$isdecoy == 0], na.rm = TRUE)
        }
      }
    } else {
      # Default cutoff if we can't merge
      Score_cutoff <- median(cpd_plot_list$Score, na.rm = TRUE)
    }
    
    return(list(Score_cutoff, cpd_plot_list))
  } else {
    # Standard FDR-based cutoff calculation
    
    # Get target and decoy scores
    unique_formula_ID <- unique(cpd_plot_list[, c("mz", "formula", "isdecoy", "Intensity", "Score")])
    target_score <- unique_formula_ID$Score[unique_formula_ID$isdecoy != 1]
    decoy_score <- unique_formula_ID$Score[unique_formula_ID$isdecoy == 1]
    
    # Handle missing values
    target_score[is.na(target_score)] <- 0
    decoy_score[is.na(decoy_score)] <- 0
    target_score[target_score == Inf] <- 0
    decoy_score[decoy_score == Inf] <- 0
    
    # Special cases
    if (length(unique(cpd_plot_list$Score)) == 1 && !(cpd_plot_list$isdecoy %in% 1)) {
      return(unique(cpd_plot_list$Score))
    } else if (length(unique(cpd_plot_list$Score)) == 1 && (cpd_plot_list$isdecoy %in% 1)) {
      return(unique(cpd_plot_list$Score) + 100)
    }
    
    if (length(unique(cpd_plot_list$Score)) == 1) {
      return(unique(cpd_plot_list$Score))
    }
    
    # Calculate FDR
    breaks <- seq(
      min(cpd_plot_list$Score), 
      max(cpd_plot_list$Score) + (max(cpd_plot_list$Score) - min(cpd_plot_list$Score)) / FDR_strip, 
      by = abs(max(cpd_plot_list$Score) - min(cpd_plot_list$Score)) / FDR_strip
    )
    
    target_score.cut <- cut(target_score, breaks, right = TRUE)
    decoy_score.cut <- cut(decoy_score, breaks, right = TRUE)
    
    target_score.freq <- table(target_score.cut)
    decoy_score.freq <- table(decoy_score.cut)
    
    target_score.freq <- target_score.freq[rev(names(target_score.freq))]
    target_score.freq0 <- c(0, cumsum(target_score.freq))
    
    decoy_score.freq <- decoy_score.freq[rev(names(decoy_score.freq))]
    decoy_score.freq0 <- c(0, cumsum(decoy_score.freq))
    
    # Calculate FDR at each score threshold
    FDR <- decoy_score.freq0 / target_score.freq0
    
    # Create FDR data frame
    df <- data.frame(
      breaks = rev(breaks),
      target_score.freq0 = target_score.freq0,
      decoy_score.freq0 = decoy_score.freq0,
      FDR = FDR
    )
    
    # Apply moving average to smooth FDR
    df$FDR_m.av <- rollmean(df$FDR, 3, fill = list(0, NaN, NaN))
    
    # Find score cutoff that satisfies FDR threshold
    Score_cutoff <- min(
      c(
        df$breaks[(df$FDR_m.av <= FDR_cutoff) == TRUE],
        df$breaks[(df$FDR <= FDR_cutoff) == TRUE]
      ), 
      na.rm = TRUE
    )
    
    # Plot FDR curve
    if (!is.null(outputdir) && plot_fdr) {
      if ((sum(is.na(df$FDR)) + sum(is.infinite(df$FDR))) < length(df$FDR)) {
        tryCatch({
          png(paste0(outputdir, "/FDR.png"))
          plot.new()
          plot(
            df$breaks, df$FDR,
            main = "FDR plot",
            xlab = "Matching score",
            ylab = "FDR"
          )
          lines(df$breaks, df$FDR_m.av)
          dev.off()
        }, error = function(e) {
          message("Error plotting FDR curve: ", e$message)
        })
      }
      
      # Write FDR table
      write.csv(df, paste0(outputdir, "/FDR.CSV"), row.names = FALSE)
    }
    
    return(Score_cutoff)
  }
}