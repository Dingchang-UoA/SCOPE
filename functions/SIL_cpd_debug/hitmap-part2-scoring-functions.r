### Functions for scoring and pattern matching

# Define function to merge isotopic patterns based on instrument resolution
isopattern_ppm_filter <- function(pattern, ppm) {
  pattern_ppm = as.numeric(as.character(pattern[, 1]))
  pattern_ppm_delta = numeric()
  
  filtered_pattern <- pattern[1, ]
  filtered_pattern <- t(data.frame(filtered_pattern, stringsAsFactors = FALSE))
  
  if (length(pattern_ppm) >= 2) {
    for (i in 1:(length(pattern_ppm) - 1)) {
      pattern_ppm_delta[i] = (pattern_ppm[i + 1] - pattern_ppm[i]) / pattern_ppm[i]
      
      if (length(pattern_ppm_delta[i] > (ppm / 1e+06)) == 1 && pattern_ppm_delta[i] > (ppm / 1000000)) {
        filtered_pattern <- rbind(filtered_pattern, pattern[i + 1, ])
      } else {
        newline = c((filtered_pattern[nrow(filtered_pattern), 1] * filtered_pattern[nrow(filtered_pattern), 2] + 
                    pattern[i + 1, 1] * pattern[i + 1, 2]) / 
                    (sum(filtered_pattern[nrow(filtered_pattern), 2], pattern[i + 1, 2])),
                    filtered_pattern[nrow(filtered_pattern), 2] + pattern[i + 1, 2])
        filtered_pattern[nrow(filtered_pattern), ] <- newline
      }
    }
  }
  
  matrix(filtered_pattern[, 1:2], ncol = 2) -> filtered_pattern
  rownames(filtered_pattern) = 1:nrow(filtered_pattern)
  
  return(filtered_pattern)
}

# Define function to filter peaklist
isopattern_ppm_filter_peaklist <- function(pattern, ppm, threshold = 0.001, verbose = FALSE) {
  org_feature = nrow(pattern)
  pattern <- as.data.frame(pattern)
  pattern_ppm = as.numeric(as.character(pattern[, 1]))
  pattern_ppm_delta = numeric()
  filtered_pattern <- pattern[1, ]
  
  for (i in 1:(length(pattern_ppm) - 1)) {
    pattern_ppm_delta[i] = (pattern_ppm[i + 1] - pattern_ppm[i]) / pattern_ppm[i]
    
    if (pattern_ppm_delta[i] > (ppm / 1000000)) {
      filtered_pattern <- rbind(filtered_pattern, pattern[i + 1, ])
    } else {
      newline = c((filtered_pattern[nrow(filtered_pattern), 1] * filtered_pattern[nrow(filtered_pattern), 2] + 
                  pattern[i + 1, 1] * pattern[i + 1, 2]) / 
                  (sum(filtered_pattern[nrow(filtered_pattern), 2], pattern[i + 1, 2])),
                  filtered_pattern[nrow(filtered_pattern), 2] + pattern[i + 1, 2])
      filtered_pattern[nrow(filtered_pattern), ] <- newline
    }
  }
  
  filtered_pattern <- filtered_pattern[filtered_pattern[, 2] >= max(filtered_pattern[, 2]) * threshold, ]
  rownames(filtered_pattern) = 1:nrow(filtered_pattern)
  
  if (verbose == TRUE) {
    message(paste("Original Features:", org_feature, "Filtered Features:", nrow(filtered_pattern)))
  }
  
  return(filtered_pattern)
}

# Function to score spectrum matching
Spectrum_scoring <- function(spec.top, spec.bottom, t = 0.25, b = 0, top.label = NULL,
                             bottom.label = NULL, xlim = c(50, 1200), x.threshold = 0, print.alignment = FALSE,
                             print.graphic = FALSE, output.list = FALSE, score_method = "SQRT", 
                             Peak_intensity = 0, outputfile = NULL, formula = NULL, pattern_ppm = NULL, 
                             ppm_error = 0, res = 72, anno_info = "") {
  
  top_tmp <- data.frame(mz = spec.top[, 1], intensity = spec.top[, 2])
  top_tmp$normalized <- round((top_tmp$intensity / max(top_tmp$intensity)), digits = 3)
  top_plot <- data.frame(mz = top_tmp$mz, intensity = top_tmp$normalized)   # data frame for plotting spectrum
  top <- subset(top_plot, top_plot$intensity > b)   # data frame for similarity score calculation
  
  bottom_tmp <- data.frame(mz = spec.bottom[, 1], intensity = spec.bottom[, 2])
  bottom_tmp$normalized <- round((bottom_tmp$intensity / max(bottom_tmp$intensity)), digits = 3)
  bottom_plot <- data.frame(mz = bottom_tmp$mz, intensity = bottom_tmp$normalized)   # data frame for plotting spectrum
  bottom <- subset(bottom_plot, bottom_plot$intensity > b)   # data frame for similarity score calculation
  
  alignment <- merge(top, bottom, by = 1, all = TRUE)
  alignment[, c(2, 3)][is.na(alignment[, c(2, 3)])] <- 0   # convert NAs to zero (R-Help, Sept. 15, 2004, John Fox)
  names(alignment) <- c("mz", "intensity.top", "intensity.bottom")
  
  if (print.alignment == TRUE) {
    print(alignment)
  }
  
  ## similarity score calculation
  
  match_score <- sum(alignment$intensity.top > 0 & alignment$intensity.bottom > 0) / 
                sum(alignment$intensity.top > 0)
  
  if (score_method == "SQRT") {
    u <- alignment[, 2]
    v <- alignment[, 3]
    similarity_score <- -log(as.vector(sqrt(sum(((u-v))^2)) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))))
    
  } else if (score_method == "SQRTP") {
    u <- alignment[, 2]
    v <- alignment[, 3]
    similarity_score <- -log(as.vector(sqrt(sum(((u-v))^2) / (sum(u^2) * sum(v^2)))))
    match_score <- log(sum(alignment$intensity.top > 0 & alignment$intensity.bottom > 0) / 
                    sum(alignment$intensity.top > 0))
    similarity_score <- similarity_score + match_score
    
  } else if (score_method == "SQRTC") {
    u <- alignment[, 2]
    v <- alignment[, 3]
    similarity_score <- -log(as.vector(sqrt(sum(((u-v))^2) / (sum(u^2) * sum(v^2)))))
    match_score <- (sum(alignment$intensity.top > 0 & alignment$intensity.bottom > 0) / 
                 sum(alignment$intensity.top > 0))
    similarity_score <- similarity_score * match_score
    
  } else if (score_method == "balanced-SQRT") {
    u <- alignment[, 2]
    v <- alignment[, 3]
    if (Peak_intensity < 10) (Peak_intensity = 10)
    
    similarity_score <- -log(as.vector(sqrt(sum(((u-v)*u)^2)) / (sqrt(sum(u^2)) * sqrt(sum(v^2))))) * 
                       log(Peak_intensity, 10)
    
  } else if (score_method == "Equal-SQRT") {
    u <- alignment[, 2]
    v <- alignment[, 3]
    score = -log(v / u)
    score = score[score != Inf]
    similarity_score = -log(sum(abs(score)) / length(score))
    
  } else if (score_method == "Equal-intensity-SQRT") {
    u <- alignment[, 2]
    v <- alignment[, 3]
    similarity_score <- as.vector((u %*% v) * length(v) / 
                                (sqrt(sum(u^2)) * sqrt(sum(v^2)) * log(Peak_intensity)))
  } else if (score_method == "Mix-SQRT") {
    alignment_0 = ifelse(alignment == 0, 0, 1)
    matching_score = sum(alignment_0[, 2] == alignment_0[, 3]) / nrow(alignment_0)
    alignment1 <- alignment[alignment_0[, 2] == alignment_0[, 3], ]
    
    u <- alignment1[, 2]
    v <- alignment1[, 3]
    
    similarity_score <- as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2))))
  }
  
  ## Generate plot
  similarity_score = round((similarity_score - ppm_error), digits = 7)
  if (print.graphic == TRUE) {
    suppressMessages(suppressWarnings(require(ggplot2)))
    tempfilename = outputfile
    if (!is.null(dev.list())) dev.off()
    png(tempfilename, bg = "white", width = 6.67, height = 6.67, res = res, units = "in")
    
    plot.new()
    plot.window(xlim = xlim, ylim = c(-125, 125))
    ticks <- c(-100, -50, 0, 50, 100)
    for (i in 1:length(top_plot$mz)) {
      lines(rep(top_plot$mz[i], 2), c(0, top_plot$intensity[i] * 100), col = "blue")
      points(top_plot$mz[i], top_plot$intensity[i] * 100, col = "darkblue")
    }
    for (i in 1:length(bottom_plot$mz)) {
      lines(rep(bottom_plot$mz[i], 2), c(0, -bottom_plot$intensity[i] * 100), col = "red")
      if (bottom_plot$intensity[i] != 0) {
        points(bottom_plot$mz[i], -bottom_plot$intensity[i] * 100, col = "darkred")
        if (length(pattern_ppm$plotppm[pattern_ppm$mz == bottom_plot$mz[i]]) == 1) {
          text(bottom_plot$mz[i], -bottom_plot$intensity[i] * 100 - 7, 
               pattern_ppm$plotppm[pattern_ppm$mz == bottom_plot$mz[i]], cex = 0.6, col = "gray21")
        }
      }
    }
    
    axis(2, at = ticks, labels = abs(ticks), pos = xlim[1], ylab = "Intensity")
    axis(1, pos = -125)
    lines(xlim, c(0, 0))
    rect(xlim[1], -125, xlim[2], 125)
    mtext("m/z", side = 1, line = 2)
    mtext("intensity (%)", side = 2, line = 2)
    plot.window(xlim = c(0, 20), ylim = c(-10, 10))
    text(10, 9, top.label)
    text(10, -9, bottom.label)
    
    formula_label = formula
    
    if (anno_info == "") {
      mtext(paste("Formula:", formula_label, "Score:", round(similarity_score, digits = 2), 
                  "Method:", score_method), cex = 0.9)
    } else {
      mtext(paste(anno_info, "\nFormula:", formula_label, "Score:", 
                  round(similarity_score, digits = 2), "Method:", score_method), cex = 0.9)
    }
    
    if (!is.null(pattern_ppm$norm_ppm)) {
      text(pattern_ppm[, 1], rep(-8.5, nrow(pattern_ppm)), pattern_ppm$norm_ppm)
    }
    
    dev.off()
  }
  
  if (sum(abs(similarity_score) == Inf, is.nan(similarity_score), is.na(similarity_score), na.rm = TRUE)) 
    similarity_score = 0
  
  if (output.list == FALSE) {
    return(similarity_score)
  }
  
  if (output.list == TRUE) {
    return(list(similarity.score = similarity_score, plot = p))
  }
}

# Main scoring function for compounds
SCORE_CPD <- function(formula, peaklist, isotopes = NULL, threshold = 0.001, charge = 1, ppm = 5,
                      print.graphic = FALSE, output.list = FALSE, outputfile = NULL, 
                      score_method = "SQRTP", similarity_plot_res = 72, anno_info = "", 
                      isolabel = "Normal", output_monomz = FALSE) {
  suppressMessages(suppressWarnings(require(rcdk)))
  suppressMessages(suppressWarnings(require(rcdklibs)))
  suppressMessages(suppressWarnings(require(OrgMassSpecR)))
  suppressMessages(suppressWarnings(require(enviPat)))
  suppressMessages(suppressWarnings(require(data.table)))
  suppressMessages(suppressWarnings(if (!require(rJava)) {
    if (Sys.getenv("JAVA_HOME") != "")
      Sys.setenv(JAVA_HOME = "")
    library(rJava)
  }))
  suppressMessages(suppressWarnings(require(grid)))
  
  if (is.null(isotopes)) {
    data("isotopes")
  }
  
  formula <- as.character(formula)
  
  # Generate the isotopic pattern
  pattern <- isopattern(
    isotopes,
    formula,
    threshold = threshold,
    plotit = FALSE,
    charge = charge,
    emass = 0.00054858,
    algo = 1,
    verbose = FALSE
  )
  
  if (ppm >= 25) {
    instrument_ppm = 50
  } else {
    instrument_ppm = max(c(8, ppm), na.rm = TRUE)
  }
  
  # Filter and merge the isotopic pattern using instrument resolution
  pattern = pattern[[formula]]
  pattern = isopattern_ppm_filter(pattern = matrix(pattern[, 1:2], ncol = 2), ppm = instrument_ppm)
  pattern = matrix(pattern[topN_feature(pattern[, 2], 10), ], ncol = 2)
  pattern <- matrix(pattern[, 1:2], ncol = 2)
  
  # Calculate spectrum intensity
  spectrumintensity = unlist(lapply(1:nrow(pattern), function(x, pattern, peaklist, ppm) {
    PMF_spectrum_intensity <- as.numeric(peaklist[
      peaklist$m.z >= pattern[x, 1] * (1 - ppm / 1000000) & 
      peaklist$m.z <= pattern[x, 1] * (1 + ppm / 1000000), 2])
    if (sum(!is.na(PMF_spectrum_intensity), length(PMF_spectrum_intensity) != 0) > 1) 
      sum(PMF_spectrum_intensity, na.rm = TRUE) 
    else 
      0
  }, pattern, peaklist, instrument_ppm))
  
  spectrum = data.frame(mz = pattern[, 1], Intensity = spectrumintensity)
  Peak_intensity <- max(spectrum$Intensity)
  
  # Get spectrum peaks
  spectrum_pk <- lapply(1:nrow(pattern), function(x, pattern, peaklist, instrument_ppm) {
    peaklist[peaklist$m.z >= pattern[x, 1] * (1 - (instrument_ppm) / 1000000) & 
            peaklist$m.z <= pattern[x, 1] * (1 + (instrument_ppm) / 1000000), ]
  }, pattern, peaklist, instrument_ppm)
  
  spectrum_pk <- do.call(rbind, spectrum_pk)
  spectrum_pk <- unique(spectrum_pk)
  
  if (nrow(spectrum_pk) > 1) {
    spectrum_pk = isopattern_ppm_filter_peaklist(pattern = spectrum_pk, ppm = ppm)
  }
  
  # Calculate ppm error
  pattern_ppm <- do.call(rbind, (lapply(1:nrow(pattern), function(x, pattern, spectrum_pk, instrument_ppm) {
    PMF_spectrum <- spectrum_pk[
      spectrum_pk$m.z >= pattern[x, 1] * (1 - instrument_ppm / 1000000) & 
      spectrum_pk$m.z <= pattern[x, 1] * (1 + instrument_ppm / 1000000), ]
    
    if (nrow(PMF_spectrum) >= 1) {
      spectrummz <- sum(PMF_spectrum[, 1] * PMF_spectrum[, 2]) / sum(PMF_spectrum[, 2])
      return(data.frame(mz = pattern[x, 1], delta_ppm = (spectrummz - pattern[x, 1]) / pattern[x, 1] * 1000000))
    }
    return(NULL)
  }, pattern, spectrum_pk, instrument_ppm)))
  
  if (!is.null(pattern_ppm) && nrow(spectrum_pk) > 1 && nrow(pattern_ppm) > 1) {
    pattern_ppm$norm_ppm <- pattern_ppm$delta_ppm - mean(pattern_ppm$delta_ppm)
    meanppm = mean(abs(pattern_ppm$norm_ppm))
    pattern_ppm$plotppm <- ifelse(abs(pattern_ppm$norm_ppm) > ppm, 
                               paste0(">", ppm), 
                               round((pattern_ppm$norm_ppm), digits = 1))
  } else if (!is.null(pattern_ppm) && nrow(spectrum_pk) == 1 && nrow(pattern_ppm) == 1) {
    pattern_ppm$norm_ppm <- pattern_ppm$delta_ppm
    meanppm = mean(abs(pattern_ppm$norm_ppm))
    pattern_ppm$plotppm <- ifelse(abs(pattern_ppm$norm_ppm) > ppm, 
                               paste0(">", ppm), 
                               round((pattern_ppm$norm_ppm), digits = 1))
  } else {
    meanppm = instrument_ppm
  }
  
  ppm_error = abs(pnorm(meanppm / ppm) - 0.5)
  
  # Calculate final score
  finalscore = Spectrum_scoring(
    spec.top = pattern,
    spec.bottom = spectrum,
    b = 0,
    t = mean(pattern[, 1]) * ppm / 1000000,
    top.label = "Theoretical",
    score_method = score_method,
    bottom.label = "Observed spectrum",
    xlim = c(min(range(pattern[, 1], na.rm = TRUE)) - 1, max(range(pattern[, 1], na.rm = TRUE)) + 1),
    output.list = output.list,
    print.graphic = print.graphic,
    outputfile = outputfile,
    Peak_intensity = Peak_intensity,
    res = similarity_plot_res,
    formula = formula,
    pattern_ppm = pattern_ppm,
    ppm_error = ppm_error,
    anno_info = anno_info
  )
  
  Peak_intensity <- sum(spectrum$Intensity)
  return(as.numeric(data.frame(finalscore, meanppm, Peak_intensity)))
}
