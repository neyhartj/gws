#' Measure the minor allele frequency of a SNP
#'
#'
measure.maf <- function(x, type = "numeric") {

  if (type == "numeric") {
    snp.recode <- as.numeric(x) + 1
    # Remove NA
    snp.recode <- na.omit(snp.recode)
    freq <- sum(snp.recode) / (2 * length(snp.recode))
    return(min(freq, 1 - freq))
  }

  if (type == "nucleotide") {
    snp.recode <- x[x != "NN"]
    snp.recode <- unlist(strsplit(snp.recode, split = ""))
    allele.table <- table(snp.recode)
    freq <- allele.table / sum(allele.table)
    return(min(freq, 1 - freq))
  }
}
