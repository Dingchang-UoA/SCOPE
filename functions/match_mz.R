match_mz  <- function(x,
                    ref,
                    tol = 1e-10,
                    unit = c("Da", "ppm")) {
  sref = sort(ref)
  
  i = findInterval(x, sref, all.inside = TRUE)
  
  dif1 = abs(x - sref[i])
  
  dif2 = abs(x - sref[i + 1])
  
  dif = dif1 > dif2
  
  dif1[dif] = dif2[dif]
  
  i[dif] = i[dif] + 1
  
  res_indx <- rep(NA, length(x))
  
  if (unit[1] == "ppm") {
    tol_ppm <- tol
    tol <- tol_ppm * x / 1000000
  }
  res_indx[dif1 <= tol] <- sref[i[dif1 <= tol]]
  res_indx
}