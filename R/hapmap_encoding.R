#' Change the encoding of a hapmap file
#'
#' @description
#' Changes the encoding of a TASSEL hapmap to a rrBLUP hapmap and vice-versa
#'
#' @param hapmap A TASSEL-encoded or rrBLUP-encoded hapmap file read in to R as
#' a data.frame. See \code{Details} for information on file format.
#' @param encoding The desired output encoding. Either \code{"rrBLUP"} or
#' \code{"TASSEL"}. See \code{Details} for information on file format
#'
#' @details
#' The TASSEL format is as such:
#' The first row is column names. The first 4 columns are marker name, alleles,
#' chromosome, and position, respectively. The next 7 column are additional information for
#' TASSEL. The remaining columns are samples. Genotypes are encoded in
#' diploid format (i.e. AA, AC, CC) with "NN" denoting missing data.
#'
#' The rrBLUP format is as such:
#' The first row is column names. The first 4 columns are marker name, alleles,
#' chromosome, and position, respectively. The next 7 column are additional information for
#' TASSEL. The remaining columns are samples. Genotypes are encoded in {1, 0, -1}
#' format where 1 is homozygous for the first allele, 0 is heterozygous, and -1 is
#' homozygous for the second allele. Missing data is denoted with \code{NA}.
#'
#' @return A \code{data.frame} of a hapmap encoded in the designated format.
#'
#' @export
#'
#'
hapmap_encoding <- function(hapmap, encoding = "rrBLUP") {

  # Error
  hapmap <- as.data.frame(hapmap)

  if (!encoding %in% c("TASSEL", "rrBLUP")) stop("Encoding must be 'TASSEL' or 'rrBLUP.'")

  # Check which encoding the hapmap is in
  if (class(hapmap[,12]) == "integer") {
    hapmap.class = "rrBLUP"
  }
  if (class(hapmap[,12]) == "character") {
    hapmap.class = "TASSEL"

    # Verify the hapmap is in diploid format
    if (nchar(hapmap[1,12]) == 1) stop ("The TASSEL hapmap is not in diploid format.")

  }

  if (hapmap.class == encoding) stop("The requested encoding already matches the class of the provided hapmap.")

  ### Change the encoding
  if (encoding == "TASSEL") {

    # Apply over rows
    genotype.recode <- t(apply(X = hapmap, MARGIN = 1, FUN = function(site) {

      # Grab the alleles and split
      alleles <- unlist(strsplit(x = site[2], split = "/"))

      genotype.calls <- site[-c(1:11)]

      # Iterate over genotypes
      unlist(sapply(X = genotype.calls, FUN = function(call) {
        if (is.na(call)) return("NN")
        if (call == " 1") return(paste(alleles[1], alleles[1], sep = ""))
        if (call == "-1") return(paste(alleles[2], alleles[2], sep = ""))
        if (call == " 0") return(paste(alleles[1], alleles[2], sep = ""))
      })) }))

    # Combine data and return
    return(cbind( hapmap[,c(1:11)], as.data.frame(genotype.recode) ))

  }

  if (encoding == "rrBLUP") {

    # Apply over rows
    genotype.recode <- t(apply(X = hapmap, MARGIN = 1, FUN = function(site) {

      alleles <- unlist(strsplit(x = as.character(site[2]), split = "/"))

      genotype.calls <- site[-c(1:11)]

      # Iterate over genotypes
      unlist(sapply(X = genotype.calls, FUN = function(call) {
        if (call == "NN") return(NA)
        else {

          # Split up the call
          call.split <- unlist(strsplit(as.character(call), ""))
          # Calls are converted to the first allele (i.e. reference) such that
          ## homozygous reference is 1, het is 0, homozygous alt is -1
          if (all(call.split == c(alleles[1], alleles[1]))) return(1)
          if (all(call.split == c(alleles[1], alleles[2]))) return(0)
          if (all(call.split == c(alleles[2], alleles[1]))) return(0)
          if (all(call.split == c(alleles[2], alleles[2]))) return(-1)
        }
      })) }))

    # Combine the matrix
    return(cbind( hapmap[,c(1:11)], as.data.frame(genotype.recode) ))

  }
}
