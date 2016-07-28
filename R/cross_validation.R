#' Run cross-validation of a training set
#'
#' @description
#' Perform cross-validation using a training set of genotypic and phenotypic data.
#' The function will perform cross-validation using a single vector of phenotypes (i.e.
#' from a single environment or a single trait). This procedure is useful to get
#' an idea of the accuracy of genomic prediction for different scenarios.
#'
#' @param phenotypes A n x t \code{matrix} of phenotypes where n is the number
#' of lines/entries and t is the number of traits. Row names should be entry/line
#' names and column names should be trait names. The function assumes
#' that the order of the entry/lines names is the same as in the genotype matrix.
#' No missing data allowed.
#' @param genotypes A n x m \code{matrix} of genotypic data where n is the number
#' of lines/entries and m is the number of biallelic markers. Column names should be
#' marker names and row names should be entry/line names. The matrix should be
#' coded in {1, 0, -1}, where 1 is homozygous for the first allele, 0 is heterzygous,
#' and -1 is homozygous for the second allele. The order of lines/entries should
#' match that of the phenotypes matrix. No missing data allowed.
#' @param cv.method The method of cross-validation. Choices are \code{"fractional"}
#' or \code{"k-fold"}. See \code{Details} for a description of these methods. Defaults
#' to both.
#' @param frac.train The proportion of the data to use as the training set when
#' running \code{"fractional"} cross-validation. In each repetition,
#'
#'
#'
#' @importFrom rrBLUP mixed.solve
#'
#'
#'
#'
#'
#'
cross.val <- function(phenotypes, genotypes, cv.method = c("fractional", "k-fold"), frac.train = 0.60, frac.reps = 250, folds = 10, fold.reps = 25) {

  # Error checking
  if (! cv.method %in% c("fractional", "k-fold")) stop("The 'cv.method' input is not valid. Valid inputs are 'fractional' and/or 'k-fold.'")
  if (frac.train <= 0 | frac.train >= 1) stop("'frac.train' must be between 0 and 1.")
  if (frac.reps < 1) stop("'frac.reps' must be >= 1.")
  if (folds < 2) stop("'folds' must be >= 2.")
  if (fold.reps < 1) stop("'fold.reps' must be >= 1.")

  # Data validation
  phenotypes <- as.matrix(phenotypes)
  genotypes <- as.matrix(genotypes)
  frac.train <- as.double(frac.train)
  frac.reps <- as.integer(frac.reps)
  folds <- as.integer(folds)
  fold.reps <- as.integer(fold.reps)

  # Make sure trait names are present
  traits <- colnames(phenotypes)
  if (is.null(traits)) stop("'phenotypes' must have column names.")

  # Check to make sure phenotypes and genotypes line up
  geno.lines <- row.names(genotypes)
  pheno.lines <- row.names(phenotypes)

  # Check if row names exist
  if (is.null(geno.lines)) stop("The 'genotypes' matrix requires row names.")
  if (is.null(pheno.lines)) stop("The 'phenotypes' matrix requires row names.")

  # Compare
  if (length(intersect(geno.lines, pheno.lines)) != length(geno.lines)) {

    warning("Not all of the lines in 'genotypes' match those in 'phenotypes.' Only the common lines will be used.")

    lines.to.use <- intersect(geno.lines, pheno.lines)

    # Subset
    genotypes <- as.matrix(genotypes[lines.to.use,])
    phenotypes <- as.matrix(phenotypes[lines.to.use,])
  }


  # Results list
  results <- list()

  # Perform fractional cross-validation if called
  if ("fractional" %in% cv.method) {

    # Apply a function over traits
    xval.results <- sapply(X = traits, FUN = function(trait) {
      cat("\nRunning fractional cross-validation on trait:", trait)
      frac.xval(y = as.matrix(phenotypes[,trait]), M = genotypes,
                frac.train = frac.train, reps = frac.reps) })

    # Add to the results list
    results[["fractional"]] <- xval.results

  } else { # Otherwise add NA
    results[["fractional"]] <- NA
  }

  if ("k-fold" %in% cv.method) {

    # Apply a function over traits
    xval.results <- sapply(X = traits, FUN = function(trait) {
      cat("\nRunning k-fold cross-validation on trait:", trait)
      kfold.xval(y = as.matrix(phenotypes[,trait]), M = genotypes,
                 folds = folds, reps = fold.reps) })

    # Add to the results list
    results[["k-fold"]] <- xval.results

  } else { # Otherwise add NA
    results[["k-fold"]] <- NA
  }

  # Return
  return(results)

} # Close the function
