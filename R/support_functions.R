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


#' Bootstrap a correlation
#'
#' @param data A two-column matrix of data to correlate
#' @param boot.reps An integer of how many bootstrapping replications to perform.
#' @param prob A floating point that specifies the percentile confidence interval to report from bootstrapping. Default is 0.95.
#'
#' @return A list containing the correlation coefficient (r), the standard deviation of bootstrapping correlation coefficients (r.sd.hat), and the confidence interval surrounding the mean (CI)
#'
#' @examples
#' data <- replicate(2, rnorm(100))
#' boot_cor(x = data[,1], y = data[,2], boot.reps = 1000)
#'
#' @import boot
#'
#' @export
#'
boot_cor <- function(x, y, boot.reps = 1000, alpha = 0.05) {

  # Error
  boot.reps <- as.integer(boot.reps)

  # Prob must be between 0 and 1
  alpha_check <- alpha > 0 | alpha < 1

  if (!alpha_check)
    stop("'alpha' must be between 0 and 1.")

  # Define a function for the correlation
  boot.cor <- function(input.data, i) {
    rep_data <- input.data[i,]
    return(cor(rep_data[,1], rep_data[,2]))
  }


  # First calculate the base statistic
  base_cor <- suppressWarnings(cor(x, y))

  # If the correlation is not NA, proceed
  if (!is.na(base_cor)) {

    # Perform the bootstrapping
    boot_results <- boot(data = cbind(x, y), statistic = boot.cor, R = boot.reps)

    # Standard error
    se <- sd(boot_results$t)
    # Bias
    bias <- mean(boot_results$t) - base_cor


    # Confidence interval
    ci_upper <- quantile(boot_results$t, 1 - (alpha / 2))
    ci_lower <- quantile(boot_results$t, (alpha / 2))

  } else {

    se <- bias <- ci_lower <- ci_upper <- NA

  }

  # Assemble list and return
  data.frame(cor = base_cor, se = se, bias = bias,
             ci_lower = ci_lower, ci_upper = ci_upper, row.names = NULL)
}







