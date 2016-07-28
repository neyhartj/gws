#' Fractional Cross Validation
#'
frac.xval <- function(y, M, frac.train, reps) {

  # Total number of lines
  n.lines <- nrow(M)

  # Number of training lines
  n.train <- round(frac.train * n.lines)

  # Replicate
  xval.results <- replicate(n = reps, expr = {

    # Randomly sample the index for lines to go into the training set
    train.ind <- sort(sample(x = seq(n.lines), size = n.train))
    test.ind <- setdiff(seq(n.lines), train.ind)

    # Subset the M and y inputs
    M.train <- M[train.ind,]
    M.test <- M[test.ind,]
    y.train <- y[train.ind,]
    y.test <- y[test.ind,]

    # Marker effects and GEBVs
    ans <- mixed.solve(y = y.train, Z = M.train)
    GEBV <- M.test %*% ans$u

    # Correlate
    cor(y.test, GEBV) })

  # Calculate mean and sd
  results <- list(x.bar = mean(xval.results), s = sd(xval.results))
  return(results)
} # Close the function

#' K-fold Cross Validation
#'
kfold.xval <- function(y, M, folds, reps) {

  # Total number of lines
  n.lines <- nrow(M)

  # Replicate
  xval.results <- replicate(n = reps, expr = {

    # Randomize
    lines.random <- sample(seq(n.lines))

    # Split into folds
    lines.folds <- split(x = lines.random, cut(x = seq_along(lines.random), breaks = folds))

    # Lapply over the list
    rep.results <- lapply(X = lines.folds, FUN = function(fold.ind) {

      # Designate the indices
      test.ind <- sort(fold.ind)
      train.ind <- sort(setdiff(lines.random, test.ind))

      # Subset the M and y inputs
      M.train <- M[train.ind,]
      M.test <- M[test.ind,]
      y.train <- y[train.ind,]
      y.test <- y[test.ind,]

      # Marker effects and GEBVs
      ans <- mixed.solve(y = y.train, Z = M.train)
      GEBV <- M.test %*% ans$u

      # Correlate
      cor(y.test, GEBV) })

    # Unlist
    unlist(rep.results) })

  # Vectorize and measure
  xval.results <- as.vector(xval.results)

  # Calculate mean and sd
  results <- list(x.bar = mean(xval.results), s = sd(xval.results))
  return(results)
} # Close the function







