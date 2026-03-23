### Plotting functions for HiTMaP package

# Function to plot compound matching scores
plot_cpd_matching_score <- function(cpd_plot_list, peaklist, charge = -1, ppm = 10, 
                                   outputdir = getwd(), filename_col = c("formula", "adduct"), 
                                   anno_col = filename_col, isotopes = "default", 
                                   similarity_plot_res = 72) {
  # Set instrument ppm based on input
  if (ppm >= 25) {
    instrument_ppm = 50
  } else {
    instrument_ppm = 8
  }
  
  message("Plotting matching isotopic patterns...")
  suppressMessages(suppressWarnings(require(enviPat)))
  suppressMessages(suppressWarnings(require(dplyr)))
  
  # Ensure required columns exist
  cpd_plot_list = data_test_rename(unique(c("formula", "isdecoy", "adduct", filename_col, anno_col)), cpd_plot_list)
  
  # Load isotopes data if needed
  if (identical(isotopes, "default")) {
    data("isotopes")
  }
  
  # Create decoy isotopes
  decoy_isotopes = isotopes
  decoy_isotopes[decoy_isotopes$isotope == "13C", ] = data.frame(
    element = "C", isotope = "11C", mass = 10.99664516, abundance = 0.0107, 
    ratioC = 0, stringsAsFactors = FALSE
  )
  
  # More decoy isotope modifications
  isotope_replacements <- list(
    list(from = "15N", to = list(element = "N", isotope = "13N", mass = 13.00603905, abundance = 0.00364000, ratioC = 4)),
    list(from = "2H", to = list(element = "H", isotope = "0H", mass = 0.001548286, abundance = 0.00011500, ratioC = 6)),
    list(from = "17O", to = list(element = "O", isotope = "15O", mass = 14.99069774, abundance = 0.00038000, ratioC = 3)),
    list(from = "18O", to = list(element = "O", isotope = "14O", mass = 13.99066884, abundance = 0.00205000, ratioC = 3)),
    list(from = "33S", to = list(element = "S", isotope = "31S", mass = 30.97268292, abundance = 0.0075, ratioC = 3)),
    list(from = "34S", to = list(element = "S", isotope = "30S", mass = 29.9762745, abundance = 0.0425, ratioC = 3)),
    list(from = "35S", to = list(element = "S", isotope = "29S", mass = 28.94414146, abundance = 0, ratioC = 3)),
    list(from = "36S", to = list(element = "S", isotope = "28S", mass = 27.97706058, abundance = 0.0001, ratioC = 3)),
    list(from = "37Cl", to = list(element = "Cl", isotope = "33Cl", mass = 32.965903, abundance = 0.24240000, ratioC = 3))
  )
  
  # Apply isotope replacements
  for (repl in isotope_replacements) {
    idx <- which(decoy_isotopes$isotope %in% repl$from)
    if (length(idx) > 0) {
      decoy_isotopes[idx, "element"] <- repl$to$element
      decoy_isotopes[idx, "isotope"] <- repl$to$isotope
      decoy_isotopes[idx, "mass"] <- repl$to$mass
      
      # Handle different column name options
      if ("abundance" %in% colnames(decoy_isotopes)) {
        decoy_isotopes[idx, "abundance"] <- repl$to$abundance
      } else if ("aboudance" %in% colnames(decoy_isotopes)) {
        # Handle potential typo in column name
        decoy_isotopes[idx, "aboudance"] <- repl$to$abundance
      }
      
      if ("ratioC" %in% colnames(decoy_isotopes)) {
        decoy_isotopes[idx, "ratioC"] <- repl$to$ratioC
      }
    }
  }
  
  # Create output directory if needed
  if (!dir.exists(outputdir)) {
    dir.create(outputdir, recursive = TRUE)
  }
  
  # Filter out compounds with no intensity
  cpd_plot_list <- cpd_plot_list[!is.na(cpd_plot_list$Intensity), ]
  
  # Group compounds if dplyr is available
  if (requireNamespace("dplyr", quietly = TRUE)) {
    group_cols <- unique(c("formula", "isdecoy", "adduct", filename_col, anno_col))
    # Use dynamic grouping with dplyr
    cpd_plot_list <- cpd_plot_list %>% 
      group_by_at(vars(one_of(group_cols))) %>% 
      summarise(Peptide = paste(Peptide, collapse = "_"))
  }
  
  # Process compounds
  if (nrow(cpd_plot_list) > 0) {
    cpd_plot_list_score <- vector("list", nrow(cpd_plot_list))
    
    for (i in 1:nrow(cpd_plot_list)) {
      # Prepare annotation info
      if ("Isotype" %in% anno_col) {
        non_isotype_cols <- anno_col[anno_col != "Isotype"]
        if (length(non_isotype_cols) > 0) {
          anno_info <- paste(cpd_plot_list[i, non_isotype_cols], collapse = "\n")
        } else {
          anno_info <- ""
        }
        isolabel <- as.character(cpd_plot_list[i, "Isotype"])
      } else {
        anno_info <- paste(cpd_plot_list[i, anno_col], collapse = "\n") 
        isolabel <- "Normal"
      }
      
      # Score target or decoy based on isdecoy flag
      formula_i <- as.character(cpd_plot_list$formula[i])
      
      if (cpd_plot_list$isdecoy[i] == 0) {
        # Target compound
        filename <- paste0(outputdir, "/", paste(cpd_plot_list[i, filename_col], collapse = "_"), "_target.png")
        
        cpd_plot_list_score[[i]] <- try(
          SCORE_CPD(
            formula_i,
            peaklist,
            isotopes = isotopes,
            charge = charge,
            ppm = ppm,
            isolabel = isolabel,
            print.graphic = TRUE,
            output.list = FALSE,
            outputfile = filename,
            similarity_plot_res = similarity_plot_res,
            anno_info = anno_info
          )
        )
      } else {
        # Decoy compound
        filename <- paste0(outputdir, "/", paste(cpd_plot_list[i, filename_col], collapse = "_"), "_decoy.png")
        
        cpd_plot_list_score[[i]] <- try(
          SCORE_CPD(
            formula_i,
            peaklist,
            isotopes = decoy_isotopes,
            charge = charge,
            ppm = ppm,
            isolabel = isolabel,
            print.graphic = TRUE,
            output.list = FALSE,
            outputfile = filename,
            similarity_plot_res = similarity_plot_res,
            anno_info = anno_info
          )
        )
      }
    }
    
    # Add scores to the result
    cpd_plot_list$Score <- unlist(cpd_plot_list_score)
    return(cpd_plot_list)
  } else {
    message("No compounds to plot")
    return(cpd_plot_list)
  }
}

# Function to plot compound spectrum matches
Plot_cpd_spectrum_match <- function(pimresultindex, spectrumlist, peplist, pimlist, 
                                   pimresultlist, threshold = 0.05) {
  suppressMessages(suppressWarnings(require(ggplot2)))
  suppressMessages(suppressWarnings(require(stringr)))
  
  message("Plotting significant PMF matches")
  
  # Filter results based on threshold
  pimresultsl <- pimresultindex[pimresultindex[, "Normalized Mean"] > threshold, ]
  
  # Prepare peaklist for plotting
  plotpeaklist <- spectrumlist
  plotpeaklist[, 1] <- as.numeric(plotpeaklist[, 1])
  plotpeaklist[, 2] <- as.numeric(plotpeaklist[, 2])
  plotpeaklist <- plotpeaklist[plotpeaklist[, "intensities"] > 0, ]
  
  # Create base plot
  sp <- ggplot(plotpeaklist, size = 1, aes(x = m.z, y = intensities)) + geom_line()
  
  # Process each result
  for (rownamelist in rownames(pimresultsl)) {
    # Create data frame for segments
    df <- data.frame(
      "mz" = pimlist[[rownamelist]],
      "mzend" = pimlist[[rownamelist]], 
      "yintercept" = rep(0, length(pimresultlist[[rownamelist]])), 
      "intensities" = pimresultlist[[rownamelist]],
      "peptide" = as.factor(peplist[[rownamelist]])
    )
    
    # Filter non-zero intensities
    df <- df[df[, "intensities"] > 0, ]
    
    # Create filename using helper function
    if (exists("windows_filename")) {
      filename <- windows_filename(str_split(rownamelist, "\\|")[[1]][2])
    } else {
      # Simple sanitization if helper not available
      filename <- gsub("[^a-zA-Z0-9_-]", "_", str_split(rownamelist, "\\|")[[1]][2])
    }
    
    # Create output file
    png(paste(getwd(), "/", filename, '.png', sep = ""), width = 1980, height = 1080)
    
    # Generate plot
    tempsp <- sp + 
      ggtitle(paste(rownamelist, "\nNormalized intensities:", 
                   pimresultsl[rownamelist, "Normalized Mean"])) + 
      geom_segment(
        aes(
          x = as.numeric(mz), 
          y = as.numeric(yintercept), 
          xend = as.numeric(mzend), 
          yend = as.numeric(intensities), 
          colour = peptide
        ),
        lineend = "round", 
        data = df,
        size = 1
      )
    
    print(tempsp)
    dev.off()
  }
}

# Function to remove outliers from compound scores
remove_cpd_score_outlier <- function(SMPLIST, IQR_LB = 0.75, outputdir = getwd(), abs_cutoff = -2) {
  # Check if we can bin the data
  suppressMessages(suppressWarnings(library(OneR)))
  suppressMessages(suppressWarnings(library(dplyr)))
  suppressMessages(suppressWarnings(library(ggplot2)))
  
  nbins = floor(length(SMPLIST$mz) / 500)
  
  if (nbins >= 2) {
    SMPLIST$mzbin <- bin(SMPLIST$mz, nbins = floor(length(SMPLIST$mz) / 500), method = "content")
  } else {
    message("Insufficient mz features to find the outlier")
    SMPLIST <- SMPLIST[SMPLIST$Score >= abs_cutoff, ]
    return(SMPLIST)
  }
  
  # Function to remove outliers
  remove_outliers <- function(x, na.rm = TRUE, IQR_LB = 0.75, ...) {
    qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
    H <- IQR_LB * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    # Not filtering upper outliers: y[x > (qnt[2] + H)] <- NA
    y
  }
  
  # Group by mz bin and remove outliers
  SMPLIST_bin_quantile <- SMPLIST %>% 
    group_by(mzbin) %>% 
    summarise(
      mz = mz,
      Score = remove_outliers(Score, IQR_LB = IQR_LB),
      mzbin = unique(mzbin)
    ) %>% 
    ungroup()
  
  # Remove NA scores
  SMPLIST_bin_quantile <- SMPLIST_bin_quantile[!is.na(SMPLIST_bin_quantile$Score), ]
  
  # Apply absolute cutoff
  SMPLIST_abs_cutoff <- SMPLIST_bin_quantile[SMPLIST_bin_quantile$Score < abs_cutoff, ]
  SMPLIST_bin_quantile <- SMPLIST_bin_quantile[SMPLIST_bin_quantile$Score >= abs_cutoff, ]
  
  # Create plot
  p <- ggplot(data = SMPLIST, mapping = aes(x = mz, y = Score, color = mzbin, group = mzbin)) +
    geom_point(
      data = SMPLIST,
      mapping = aes(x = mz, y = Score),
      color = "darkblue",
      group = "Cut-off",
      col = "Removed Features"
    ) +
    geom_point(
      data = SMPLIST_bin_quantile,
      mapping = aes(x = mz, y = Score),
      color = "red",
      group = "Cut-off",
      col = "Kept Features"
    ) +
    guides(color = guide_legend("m/z bin size")) +
    labs(x = 'm/z') +
    theme_classic()
  
  # Save plot if output directory is provided
  if (!is.null(outputdir)) {
    summary_dir <- paste0(outputdir, "/Summary folder")
    dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
    
    png(
      paste0(summary_dir, "/Peptide_score_outlier.png"),
      width = 1200,
      height = 800,
      res = 150
    )
    print(p)
    dev.off()
  }
  
  # Return filtered data
  return(SMPLIST[SMPLIST$mz %in% SMPLIST_bin_quantile$mz, ])
}

# Function to build adduct list - minimal implementation
Build_adduct_list <- function() {
  # Create basic adduct list for negative mode
  adduct_list <- data.frame(
    Name = c("M-H", "M+Cl", "M-H2O-H", "M+Na-2H", "M+K-2H", "2M-H", "3M-H"),
    Mass = c(-1.007825, 34.969402, -19.01839, 20.974666, 36.948606, -1.007825, -1.007825),
    Charge = c(-1, -1, -1, -1, -1, -1, -1),
    Mult = c(1, 1, 1, 1, 1, 2, 3),
    stringsAsFactors = FALSE
  )
  
  return(adduct_list)
}

# Utility function to resolve multiple modes in a distribution
resolve_multi_modes <- function(mm, adjust = 0.25) {
  suppressMessages(suppressWarnings(require(dplyr)))
  suppressMessages(suppressWarnings(require(tidyr)))
  
  # Helper function to get modes
  get.modes2 <- function(x, adjust, signifi, from, to) {
    den <- density(x, kernel = c("gaussian"), adjust = adjust, from = from, to = to)
    den.s <- smooth.spline(den$x, den$y, all.knots = TRUE, spar = 0.1)
    s.1 <- predict(den.s, den.s$x, deriv = 1)
    s.0 <- predict(den.s, den.s$x, deriv = 0)
    den.sign <- sign(s.1$y)
    a <- c(1, 1 + which(diff(den.sign) != 0))
    b <- rle(den.sign)$values
    df <- data.frame(a, b)
    df <- df[which(df$b %in% -1), ]
    modes <- s.1$x[df$a]
    density <- s.0$y[df$a]
    df2 <- data.frame(modes, density)
    df2$sig <- signif(df2$density, signifi)
    df2 <- df2[with(df2, order(-sig)), ]
    
    # Convert to tidyr format
    df <- as.data.frame(
      df2 %>%
        mutate(m = min_rank(desc(sig))) %>%
        group_by(m) %>%
        summarize(a = paste(format(round(modes, 2), nsmall = 2), collapse = ',')) %>%
        spread(m, a, sep = '')
    )
    
    colnames(df) <- paste0("m", 1:length(colnames(df)))
    return(df)
  }
  
  # Create data frame and plot
  mmdf <- data.frame(mm = mm)
  suppressMessages(suppressWarnings(require(ggplot2)))
  
  # Create density plot
  p <- ggplot(mmdf, aes(mm)) + 
    geom_density(adjust = adjust) + 
    xlim((min(mm) - 1), (max(mm) + 1))
  
  # Get modes (with significance level 2)
  modes <- get.modes2(mm, adjust = adjust, 2, min(mm) - 1, max(mm) + 1)
  
  return(list(p = p, modes = modes))
}
