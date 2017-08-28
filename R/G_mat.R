#' Calculate a genomic relationship matrix from marker information
#'
#' @description
#' Calculates the genomic relationship matrix \code{G} using methods described
#' by Van Raden (2008).
#'
#' @param X A matrix (\eqn{n \times m}) of marker data on \eqn{n} entries and
#' \eqn{m} biallelic markers, coded as {-1, 0, 1}. Fractional (imputed)
#' and missing values (NA) are allowed. Must have row (entry) names and column
#' (marker) names.
#' @param base.pop A \code{character} vector of entry names that make up the
#' base population (i.e. unselected) for calculating marker allele frequencies
#' (see \emph{Details}).
#' @param min.MAF See \code{\link[rrBLUP]{A.mat}}
#' @param max.missing See \code{\link[rrBLUP]{A.mat}}
#' @param impute.method The method for imputation. Can be "mean", "EM", or "pass".
#' See \code{\link{rrBLUP}[A.mat]} for "mean" or "EM". If "pass", markers
#' are assumed imputed.
#' @param tol See \code{\link[rrBLUP]{A.mat}}
#' @param shrink See \code{\link[rrBLUP]{A.mat}}
#'
#' @importFrom rrBLUP A.mat
#'
#' @export
#'
G_mat <- function(X, base.pop = NULL, min.MAF = NULL, max.missing = NULL,
                impute.method = c("mean", "EM", "pass"), tol = 0.02, shrink = FALSE) {

  ## Check the X input
  # It must have dimnames
  if (any(sapply(X = dimnames(X), is.null))) {
    stop("The marker matrix X must have row and column names.")
  }

  # Make sure entries in 'base.pop' are in the matrix
  if (any(!base.pop %in% row.names(X))) {
    stop("Not all entries in the input 'base.pop' are in the marker matrix X.")
  }

  impute.method <- match.arg(impute.method)

  # Impute, if not 'pass'
  if (impute.method != "pass") {

    # Use A.mat to impute
    X_impute <- A.mat(X = X, min.MAF = min.MAF, max.missing = max.missing,
                      impute.method = impute.method, tol = tol, n.core = 1,
                      shrink = shrink, return.imputed = TRUE)$imputed

  } else {
    # Make sure there are no missing elements
    if (any(is.na(X))) {
      stop("If 'impute.method' is 'pass', there can be no missing elements.")
    }

    # Re-name
    X_impute <- X

  }

  # Extract the base population
  X_0 <- X_impute[base.pop,]

  # Find the frequency of the 1 allele at each locus
  p <- colMeans(X_0 + 1) / 2

  # Calculate the P matrix
  P <- matrix(data = 2 * (p - 0.5), nrow = nrow(X_impute), ncol = ncol(X_impute), byrow = T)
  # Calculate the standardization constant
  c <- 2 * sum(p * (1 - p))

  # Subtract the P matrix from the original M matrix
  Z <- X_impute - P

  # Calculate the relationship matrix
  G <- tcrossprod(Z) / c

  # Return the matrix
  return(G)

}





