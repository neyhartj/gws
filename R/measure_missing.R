#' Measure the missingness of a SNP or entry
#'
#'
measure.missing <- function(x, type = "numeric") {

  if (type == "numeric") {
    return( sum(is.na(x)) / length(x) )
  }

  if (type == "nucleotide") {
    return( sum(x == "NN") / length(x) )
  }
}
