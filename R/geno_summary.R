#' Summarize a marker matrix
#'
#' @param genos.mat A n x m marker matrix where n is the number of entries and m
#' is the number of markers.
#' @param print.plot Logical. Should plots of the statistics be displayed?
#'
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#'
#' @export
#'
geno_summary <- function(genos.mat, print.plot = TRUE) {

  # Deal with input
  genos.mat <- genos.mat %>%
    as.matrix()

  # Line names
  n.names <- rownames(genos.mat)
  # Marker names
  m.names <- colnames(genos.mat)

  ## Calculations for markers

  # Missingness
  m.missingness <- genos.mat %>%
    is.na() %>%
    colMeans()

  # MAF
  m.maf <- (genos.mat + 1) %>%
    colMeans(na.rm = T) %>%
    {. / 2} %>%
    sapply(FUN = function(freq) min(freq, 1 - freq))

  # Heterozygosity
  m.het <- genos.mat %>%
    {. == 0} %>%
    colMeans(na.rm = T)

  ## Calculations for entries

  # Missingness
  n.missingness <- genos.mat %>%
    is.na() %>%
    rowMeans()

  # Heterozygosity
  n.het <- genos.mat %>%
    {. == 0} %>%
    rowMeans(na.rm = T)

  ## Assemble data.frame
  m.df <- data.frame(
    marker = m.names,
    p.missing = m.missingness,
    maf = m.maf,
    p.het = m.het
  ) %>%
    tbl_df() %>%
    gather(metric, value, -marker)

  n.df <- data.frame(
    entry = n.names,
    p.missing = n.missingness,
    p.het = n.het
  ) %>%
    tbl_df() %>%
    gather(metric, value, -entry)


  ## Plot, if desired
  if (print.plot) {

    gp <- ggplot(data = m.df, aes(x = value)) +
      geom_histogram(bins = 15) +
      facet_grid(~ metric, scales = "free_x")

    print(gp)

  }

  # Return the information
  return(list(
    m.df = m.df,
    n.df = n.df
  ))

} # Close the function




