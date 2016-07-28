#' Measure the heterozygosity of a SNP or entry
#'
#'
measure.het <- function(x, type = "numeric") {

  if (type == "numeric") {
    snp.recode <- na.omit(x)
    return( sum(snp.recode == 0) / length(snp.recode) )
  }

  if (type == "nucleotide") {
    snp.recode <- x[x != "NN"]
    n.het <- sum(sapply(X = strsplit(snp.recode, ""), FUN = function(site) length(unique(site)) == 2))
    return( n.het / length(snp.recode) )
  }
}
