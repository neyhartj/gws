#' Filter a hapmap or genotype matrix
#'
#' @description
#' Filters a hapmap or genotype matrix based on user-defined limits of SNP missingness,
#' SNP minor allele frequency, SNP heterozygosity, entry missingness, and entry
#' heterozygosity.
#'
#' @param x A hapmap or genotype matrix. Heuristics are used to determine the format.
#' Both may be encoded using TASSEL or rrBLUP standards. See \code{Details} for
#' information on file format.
#' @param min.maf The minimum minor allele frequency cutoff to keep a SNP.
#' @param min.snp.missing The maxmimum missingness proportion
#'
#'
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
#' @import dplyr
#' @import stringr
#' @import tidyr
#'
#' @export
#'
#'
filter_genos <- function(genos, min.maf = 0, max.mar.missing = 1,
                         max.entry.missing = 1, max.mar.het = 1, max.entry.het = 1,
                         print.plot = FALSE, verbose = TRUE) {

  # Separate the input
  genos.components <- convert_hapmap(genos)

  # Separate
  genos.mat <- genos.components$marker.genos
  genos.info <- genos.components$geno.info

  # Count markers and lines
  n <- nrow(genos.mat) # Number of lines
  m <- ncol(genos.mat) # Number of markers

  # Print notification
  if (verbose)
    cat(str_c("\nThere are ", n, " lines and ", m, " markers in the input genotype data."))

  # create results list
  results.out <- list()

  ## Start calculations
  geno.stats <- geno_summary(genos.mat = genos.mat, print.plot = F)

  # Separate and expand
  m.df <- geno.stats$m.df %>%
    spread(metric, value)

  # Filter markers first
  m.df1 <- m.df %>%
    filter(maf >= min.maf & p.het <= max.mar.het & p.missing <= max.mar.missing)

  # Filtered marker names
  m.new <- m.df1 %>%
    select(marker) %>%
    unlist()

  # subset the geno matrix
  genos.mat1 <- genos.mat %>%
    subset.matrix(select = m.new)

  # Re-run the stats
  geno.stats1 <- geno_summary(genos.mat = genos.mat1, print.plot = F)

  # Separate and expand
  n.df <- geno.stats1$n.df %>%
    spread(metric, value)

  # Filter entries
  n.df1 <- n.df %>%
    filter(p.het <= max.entry.het & p.missing <= max.entry.missing)

  # Filtered marker names
  n.new <- n.df1 %>%
    select(entry) %>%
    unlist()

  # Refilter the marker matrix
  genos.mat2 <- genos.mat1 %>%
    subset.matrix(subset = row.names(.) %in% n.new)

  # Re-run stats
  geno.stats2 <- geno_summary(genos.mat = genos.mat2, print.plot = F)



  ## Remove markers first
  # Compute marker stats
  m.stats <- marker.stats(geno.in = geno.filt)
  # Decide which markers to keep
  keep.m.MAF <- names(which(m.stats$MAF >= min.maf))
  keep.m.het <- names(which(m.stats$marker.het <= max.mar.het))
  keep.m.miss <- names(which(m.stats$marker.miss <= max.mar.missing))
  keep.m.list <- list(keep.m.MAF, keep.m.het, keep.m.miss)
  # Make a master index of those markers that pass all filters
  m.to.keep <- Reduce(f = intersect, x = keep.m.list)
  # Filter
  geno.filt <- geno.filt[,m.to.keep]

  ## Remove entries second
  # Calculate entry stats
  e.stats <- entry.stats(geno.in = geno.filt)
  # Decide which entries to keep
  keep.e.miss <- which(e.stats$entry.miss <= max.entry.missing)
  keep.e.het <- which(e.stats$entry.het <= max.entry.het)
  keep.e.list <- list(keep.e.miss, keep.e.het)
  e.to.keep <- Reduce(f = intersect, x = keep.e.list)
  # Filter
  geno.filt <- geno.filt[e.to.keep,]

  # Notify
  if (verbose) {
    cat( paste( "\n\n", (m - length(keep.m.MAF)), " markers did not pass the minor allele frequency filter.", sep = "" ) )
    cat( paste( "\n", (m - length(keep.m.het)), " markers did not pass the marker heterozygosity filter.", sep = "" ) )
    cat( paste( "\n", (m - length(keep.m.miss)), " markers did not pass the marker missingness filter.", sep = "" ) )
    cat(paste("\n",ncol(geno.filt), " markers were retained from an original ", m, " markers (", round(100*(ncol(geno.filt)/m), 3), "%).", sep = ""))

    cat( paste( "\n\n", (n - length(keep.e.het)), " entries did not pass the entry heterozygosity filter.", sep = "" ) )
    cat( paste( "\n", (n - length(keep.e.miss)), " entries did not pass the entry missingness filter.", sep = "" ) )
    cat(paste("\n",nrow(geno.filt), " entries were retained from an original ", n, " lines (", round(100*(nrow(geno.filt)/n), 3), "%).", sep = ""))
  }

  # New dimension calculation
  n <- nrow(geno.filt)
  m <- ncol(geno.filt)


  # Put the genotype data into the list
  results.out[[1]] <- geno.filt; names(results.out)[1] <- "geno.out"

  # Create a genotype table to output
  g.table <- table(geno.filt)/(nrow(geno.filt)*ncol(geno.filt))
  results.out[[2]] <- g.table; names(results.out)[2] <- "g.table"

  # Output graphs if told to
  if (graphs) {
    # Recalculate marker parameters post-filtering
    m.stats <- marker.stats(geno.in = geno.filt)
    e.stats <- entry.stats(geno.in = geno.filt)

    # Plot
    par(mfrow = c(2,3))
    # Site-frequency spectrum
    hist(m.stats$MAF, main = "Site Frequency Spectrum", xlab = "Minor Allele Frequency", xlim = c(0,0.5))
    # Missing data per marker
    hist(m.stats$marker.miss, main = "Frequency of Missing Marker Data", xlab = "Marker Missing Data Proportion", xlim = c(0,max.mar.missing))
    # Marker heterozygosity
    hist(m.stats$marker.het, main = "Frequency of Marker Heterozygosity", xlab = "Marker Heterozygosity", xlim = c(0,max.mar.het))
    # Blank plot to take up space
    plot.new()
    # Missing data per line
    hist(e.stats$entry.miss, main = "Frequency of Missing Entry Data", xlab = "Line Missing Data Proportion", xlim = c(0,max.entry.missing))
    # Line heterozygosity
    hist(e.stats$entry.het, main = "Frequency of Entry Heterozygosity", xlab = "Entry Heterozygosity Proportion", xlim = c(0,max.entry.het))
  }

  # Return the variable
  return(results.out)
}
