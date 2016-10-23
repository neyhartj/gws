#' Convert a hapmap to a genotype matrix
#'
#' @param hapmap A hapmap data.frame using the rrBLUP or TASSEL encoding.
#'
#' @import dplyr
#' @import stringr
#'
convert_hapmap <- function(hapmap) {

  # Force to tibble
  hapmap1 <- hapmap %>%
    tbl_df()

  # Determine if rrBLUP or TASSEL
  rrblup.cols <- "rs|rs#|alleles|chrom|pos"
  tassel.cols <- "rs|rs#|alleles|chrom|pos|assembly#|protLSID|assayLSID|panel|QCcode"

  rrblup.detect <- genos1 %>%
    names() %>%
    str_detect(pattern = "rs|rs#|alleles|chrom|pos") %>%
    sum()

  tassel.detect <- genos1 %>%
    names() %>%
    str_detect(pattern = tassel.cols) %>%
    sum()

  # Determine the number of matches
  if (rrblup.detect != 4 & tassel.detect != 11) {
    type = "geno_mat"

  } else {

    if (rrblup.detect == 4) {
      type = "rrblup"

    }
    if (tassel.detect == 11) {
      type = "tassel"

    }
  }

  # Return the matrix if the type is geno_mat
  if (type == "geno_mat")
    return(list(geno.info = NA,
                marker.genos = as.matrix(genos1)))

  # Separate based on type
  if (type == "rrblup") {

     # Separate the snp info
     snp.info <- genos1 %>%
       select(1:4)

     # Transpose the marker matrix
     marker.genos <- genos1 %>%
       select(-1:-4) %>%
       t()

     colnames(marker.genos) <- genos1$rs

  }

  if (type == "tassel") {

    # Separate the snp info
    snp.info <- genos1 %>%
      select(1:11)

    # Transpose the marker matrix
    marker.genos <- genos1 %>%
      select(-1:-11) %>%
      t()

    colnames(marker.genos) <- genos1$rs

  }

  # Build and return a list
  return(list(geno.info = snp.info,
              marker.genos = marker.genos))

} # Close the function






