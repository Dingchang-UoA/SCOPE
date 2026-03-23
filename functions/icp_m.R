# This function is modified based on LOMAR::icp(), which can compute the transformed points within the 
# given geometry. X is the model geometry, Y is the observed geometry and Y_full is the combination 
# of outline and inner points of the observed geometry 

icp_m <- function (X, Y, Y_full =NULL, weights = NULL, iterations = 1000, subsample = NULL, 
                   scale = TRUE, tol = 0.0001) 
{
  if (dim(Y)[2] != dim(X)[2]) {
    stop("ERROR: Point sets must have same dimension")
  }
  if (!is.null(subsample) && (subsample >= 1 || subsample <= 
                              0)) {
    stop("ERROR: subsample must be a fraction between 0 and 1")
  }
  D <- dim(X)[2]
  X <- as.matrix(X)
  Xs <-standardize_coordinates(X)
  X <- Xs[["X"]]
  Y <- as.matrix(Y)
  Ys <- standardize_coordinates(Y)
  Y <- Ys[["X"]]
  
  Y_full <-as.matrix(Y_full)
  Y_full_s <- t(apply(Y_full,1, function(x){
    x - Ys$mu
  }))
  Y_full_s <- Y_full_s/Ys$sigma
  Y_full_t <- Y_full_s
  
  Y0 <- NULL
  if (!is.null(subsample)) {
    sample.X <- sample(nrow(X), floor(subsample * nrow(X)))
    sample.Y <- sample(nrow(Y), floor(subsample * nrow(Y)))
    X <- X[sample.X, ]
    Y0 <- Y
    Y <- Y[sample.Y, ]
    if (!is.null(weights)) {
      weights <- weights[sample.Y, sample.X]
    }
  }
  M <- dim(Y)[1]
  N <- dim(X)[1]
  R <- diag(D)
  t <- rep(0, D)
  s <- 1
  Yt <- Y
  iter <- 0
  converged <- FALSE
  err <- 1e+154
  R.final <- diag(D)
  t.final <- rep(0, D)
  s.final <- 1
  
  while (iter < iterations) {
    iter <- iter + 1
    nn <- RANN::nn2(X, Yt, k = 1, searchtype = "standard")
    Xc <- X[nn[["nn.idx"]], ]
    old.err <- err
    err <- sum((Xc - Yt)^2)
    
    message(sprintf("Iteration: %d, Distance: %f", iter, err))
    
    if (abs(old.err - err)/err < tol) {
      converged <- TRUE
      break
    }
    muX <- colMeans(Xc)
    muY <- colMeans(Yt)
    Xhat <- as.matrix(Xc - muX)
    Yhat <- as.matrix(Yt - muY)
    A <- t(Xhat) %*% Yhat
    svd.A <- svd(A)
    C <- diag(D)
    C[length(C)] <- det(svd.A$u) * det(svd.A$v)
    R <- svd.A$u %*% C %*% t(svd.A$v)
    if (scale) {
      s <- sum(A * R)/sum(t(Yhat) %*% Yhat)
      s.final <- s.final * s
    }
    t <- muX - s * R %*% muY
    Yt <- apply_transformation(Yt, R, t, s)
    
    if (!is.null(Y_full)) {
      Y_full_t <- apply_transformation(Y_full_t, R, t, s)
    }
    
    R.final <- R %*% R.final
    t.final <- t.final + t
  }
  
  # 
  # if (!is.null(Y0)) {
  #   Yt <- apply_transformation(Y0, R.final, t.final, s.final)
  # }
  # Yt <- restore_coordinates(Yt, Xs[["mu"]], Ys[["sigma"]])
  # return(list(Y = Yt, R = R.final, t = t.final * Ys[["sigma"]], 
  #             s = s, iter = iter, conv = converged))
  
  return(list(X= as.data.frame(Xc) ,Y = as.data.frame(Yt), Y_full_t = Y_full_t, R = R.final, 
              scale_factor = s.final * Xs[["sigma"]] / Ys[["sigma"]],
              iter = iter, conv = converged))
}
