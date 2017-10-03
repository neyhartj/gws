#' Bootstrap a correlation
#'
#' @param data A two-column matrix of data to correlate
#' @param boot.reps An integer of how many bootstrapping replications to perform.
#' @param prob A floating point that specifies the percentile confidence interval to report from bootstrapping. Default is 0.95.
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

  # Perform the bootstrapping
  boot_results <- boot(data = cbind(x, y), statistic = boot.cor, R = boot.reps)

  # Extract the original statistic
  cor0 <- boot_results$t0
  # Standard error
  se <- sd(boot_results$t)
  # Bias
  bias <- mean(boot_results$t) - cor0
  # Confidence interval
  ci_lower <- (2 * cor0) - quantile(boot_results$t, 1 - (alpha / 2))
  ci_upper <- (2 * cor0) - quantile(boot_results$t, (alpha / 2))

  # Assemble list and return
  data.frame(cor = boot_results$t0, se = se, bias = bias,
             ci_lower = ci_lower, ci_upper = ci_upper, row.names = NULL)
}
