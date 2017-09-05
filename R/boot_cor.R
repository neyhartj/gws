#' Bootstrap a correlation
#'
#' @param data A two-column matrix of data to correlate
#' @param boot.reps An integer of how many bootstrapping replications to perform.
#' @param prob A floating point that specifies the percentile confidence interval to report from bootstrapping. Default is 0.95.
#' @return A list containing the correlation coefficient (r), the standard deviation of bootstrapping correlation coefficients (r.sd.hat), and the confidence interval surrounding the mean (CI)
#'
#' @examples
#' data <- replicate(2, rnorm(100))
#' boot.cor(data = data, boot.reps = 1000)
#'
#' @import boot
#'
#' @export
#'
boot_cor <- function(x, y, boot.reps, prob = 0.95) {

  # Error
  boot.reps <- as.integer(boot.reps)

  # Prob must be between 0 and 1
  prob_check <- prob > 0 | prob < 1

  if (!prob_check)
    stop("'prob' must be between 0 and 1.")

  # Define a function for the correlation
  boot.cor <- function(input.data, i) {
    rep_data <- input.data[i,]
    return(cor(rep_data[,1], rep_data[,2]))
  }

  # Perform the bootstrapping
  boot_results <- boot(data = cbind(x, y), statistic = boot.cor, R = boot.reps)

  # Calculate a confidence interval
  CI <- quantile(boot_results$t, probs = c(1 - prob, prob), na.rm = TRUE)

  # Assemble list and return
  data.frame(r_hat = boot_results$t0, r_sd_hat = sd(boot_results$t),
             CI_lower = CI[1], CI_upper = CI[2])
}
