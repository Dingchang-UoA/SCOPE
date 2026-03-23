### Core functions for HiTMaP package

MonoisotopicMass_sil <- function(formula = list(), isotopes = list(), charge = 0, Isotype = "Normal", c_num = 6) {
  defaultFormula <- list(C = 0, H = 0, N = 0, O = 0, S = 0, 
                         P = 0, Br = 0, Cl = 0, F = 0, Si = 0, x = 0, '[13]C' = 0)
  defaultFormula[names(formula)] <- formula
  defaultIsotopes <- list(C = 12, H = 1.0078250321, N = 14.0030740052, 
                          O = 15.9949146221, S = 31.97207069, P = 30.97376151, 
                          Br = 78.9183376, Cl = 34.96885271, F = 18.9984032, Si = 27.9769265327, 
                          x = 0, '[13]C' = 13.003355)
  
  defaultIsotopes_13c <- list(C = 13.003355, H = 1.0078250321, N = 14.0030740052, 
                          O = 15.9949146221, S = 31.97207069, P = 30.97376151, 
                          Br = 78.9183376, Cl = 34.96885271, F = 18.9984032, Si = 27.9769265327, 
                          x = 0, '[13]C' = 13.003355)
  defaultIsotopes[names(isotopes)] <- isotopes
  if(Isotype == "Normal"){
    if (charge < 0 & abs(charge) > defaultFormula$H) 
      stop("the number of negative charges exceeds the number of hydrogens in the formula list")
    
    mass <- (defaultFormula$C * defaultIsotopes$C + defaultFormula$H * 
             defaultIsotopes$H + defaultFormula$N * defaultIsotopes$N + 
             defaultFormula$O * defaultIsotopes$O + defaultFormula$S * 
             defaultIsotopes$S + defaultFormula$P * defaultIsotopes$P + 
             defaultFormula$Br * defaultIsotopes$Br + defaultFormula$Cl * 
             defaultIsotopes$Cl + defaultFormula$F * defaultIsotopes$F + 
             defaultFormula$Si * defaultIsotopes$Si + defaultFormula$x * 
             defaultIsotopes$x + defaultIsotopes$'[13]C' * defaultFormula$'[13]C')
    if (charge != 0) 
      mass <- abs((mass + charge * 1.007276466)/charge)
  } else if(Isotype == "SIL") {
    if (charge < 0 & abs(charge) > defaultFormula$H) 
      stop("the number of negative charges exceeds the number of hydrogens in the formula list")
    if (defaultFormula$C >= 6){
      mass <- (6 * defaultIsotopes_13c$C + (defaultFormula$C-6) * defaultIsotopes$C + defaultFormula$H * 
               defaultIsotopes_13c$H + defaultFormula$N * defaultIsotopes_13c$N + 
               defaultFormula$O * defaultIsotopes_13c$O + defaultFormula$S * 
               defaultIsotopes_13c$S + defaultFormula$P * defaultIsotopes_13c$P + 
               defaultFormula$Br * defaultIsotopes_13c$Br + defaultFormula$Cl * 
               defaultIsotopes_13c$Cl + defaultFormula$F * defaultIsotopes_13c$F + 
               defaultFormula$Si * defaultIsotopes_13c$Si + defaultFormula$x * 
               defaultIsotopes_13c$x)
    } else {
      mass <- (defaultFormula$C * defaultIsotopes$C + defaultFormula$H * 
               defaultIsotopes$H + defaultFormula$N * defaultIsotopes$N + 
               defaultFormula$O * defaultIsotopes$O + defaultFormula$S * 
               defaultIsotopes$S + defaultFormula$P * defaultIsotopes$P + 
               defaultFormula$Br * defaultIsotopes$Br + defaultFormula$Cl * 
               defaultIsotopes$Cl + defaultFormula$F * defaultIsotopes$F + 
               defaultFormula$Si * defaultIsotopes$Si + defaultFormula$x * 
               defaultIsotopes$x)
    }

    if (charge != 0) 
      mass <- abs((mass + charge * 1.007276466)/charge)
  }

  return(mass)
}

# Helper function to test and rename columns
data_test_rename <- function(required_col, df) {
  testcolumeresult = testcolume(df, required_col)
  
  if (length(testcolumeresult$failcol) > 0) {
    stop(paste(testcolumeresult$failcol, "column is missing in the input file, please check your datafile"))  
  }
  
  if (length(testcolumeresult$renamecol) > 0) {
    lapply(testcolumeresult$renamecol, gsub, testcolumeresult$renamecol, ignore.case = TRUE, x = names(df))
    
    for (i in length(testcolumeresult$renamecol)) {
      names(df) = gsub(testcolumeresult$renamecol, testcolumeresult$renamecol, ignore.case = TRUE, x = names(df))
    }
  }
  
  if (length(testcolumeresult$duplicatecol) > 0) {
    message(paste("Found ambiguous or duplicate column: ", testcolumeresult$duplicatecol))
  }
  
  return(df)
}

# Helper function to test columns
testcolume <- function(df, testcolnames) {
  suppressMessages(suppressWarnings(require(stringr)))
  
  testcolnames = unique(testcolnames)
  
  test = sapply(testcolnames, FUN = grepl, names(df))
  
  testresult = sapply(colnames(test), FUN = function(x, df) { sum(df[, x]) }, test) == 0
  
  case_test = sapply(testcolnames, FUN = grepl, names(df), ignore.case = TRUE)
  
  case_testresult = sapply(testcolnames, FUN = function(x, df) { sum(df[, x]) }, case_test) > 1
  
  duplicatecol = names(case_testresult)[case_testresult == TRUE]
  
  if (length(duplicatecol) == 0) { duplicatecol = NULL }
  
  testfail = names(testresult)[testresult == TRUE]
  
  passcol = names(testresult)[testresult == FALSE]
  
  renamecol = NULL
  
  if (length(passcol) < length(testcolnames)) {
    case_test = sapply(testfail, FUN = grepl, names(df), ignore.case = TRUE)
    case_testresult = sapply(colnames(case_test), FUN = function(x, df) { sum(df[, x]) }, case_test) == 1
    
    renamecol = names(case_testresult)[case_testresult == TRUE]
    
    if (length(renamecol) > 0) {
      failcol = testfail[!grep(renamecol, testfail)]  
    } else {
      failcol = NULL
    }
  } else {
    failcol = NULL
  }
  
  return(list(passcol = passcol, renamecol = renamecol, failcol = failcol, duplicatecol = duplicatecol))
}

# Helper function to convert atom lists to formula strings
atomes.to.formula <- function(unique_formula_list) {
  lapply(names(unique_formula_list), function(x) {
    unique_formula_list[[x]] -> temp
    temp[temp == 0] <- NULL
    temp_for <- paste0(names(temp), unlist(temp), collapse = "")
    return(temp_for)
  }) -> unique_formula_list_for
  return(unlist(unique_formula_list_for))
}

# Function to select top N features
topN_feature <- function(vec, n, dec = TRUE) {
  inx <- order(vec, decreasing = dec)[1:min(n, length(vec))]
  vec <- rep(FALSE, length = length(vec))
  vec[inx] <- TRUE
  return(vec)
}
