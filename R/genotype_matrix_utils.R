
# Common marker finder function
# Common marker finder function
common.markers <- function(geno.list = NULL # A list of two or more genotype matricies
) {
  # Error reporting
  if (is.null(geno.list)) {
    stop ("The list genotype matricies was not specified.")
  }

  # Extract marker names from all matricies
  m.names <- sapply(X = geno.list, FUN = function(x) return(colnames(x)) )
  # Find common markers
  m.same <- Reduce(f = intersect, x = m.names)

  cat(paste("\nThere are ", length(m.same), " markers in common.", sep = ""))

  # Subset the matricies based on those common markers
  geno.list.subset <- sapply(X = geno.list, FUN = function(mat) mat[,m.same] )

  # Name the list components
  for (i in 1:length(geno.list.subset)) {
    names(geno.list.subset)[i] <- paste("geno.mat", i, sep = "")
  }

  return(geno.list.subset)
}

## geno.matrix.stats
# Function to calculate various statistics using a genotype matrix
# Will also output graphs, if desired
geno.stats <- function(geno.in = NULL, # A hapmap file or n x m genotype matrix (see above)
                       graphs = TRUE # Logical scalar to print graphs of various statistics
                     ) {

  ## Input handling
  # Detect the type of input for g.in
  # If the column names match the expected, then the class of the input is hapmap
  if(is.character(geno.in[1,1]) && is.character(geno.in[1,2])) {
    class.in <- "hapmap"
  } else {
    class.in <- "geno.mat"
  }

  # Reformatting the g.in if it is a hapmap file
  if(class.in == "hapmap") {
    # Reformat the data geno.in a genotype matrix
    row.names(geno.in) <- geno.in[,1] # Retain the marker names in the matrix
    geno.in <- geno.in[,-c(1:4)]
    geno.in <- t(geno.in)
  }

  # Count markers and lines
  n <- nrow(geno.in) # Number of lines
  m <- ncol(geno.in) # Number of markers

  ## Stat calculator functions
  marker.stats <- function(geno.in) {
    ## Determine the minor allele frequency
    # Adding one results in each component of geno.in being the count of the number of reference alleles
    # at a SNP locus per individual
    MAF <- apply(X = geno.in + 1, MARGIN = 2, FUN = function(snp) {
      # The mean / 2 is the frequency of the reference allele
      p.ref <- mean(snp, na.rm = T) / 2
      # The MAF is the smallest of p.ref and 1 - p.ref
      MAF <- min(p.ref, (1 - p.ref))
      return(MAF)
    })
    # Find the mean
    mu.MAF <- mean(MAF, na.rm = T)

    ## Proportion of missing data
    marker.miss <- apply(X = geno.in, MARGIN = 2, FUN = function(snp) return(sum(is.na(snp)) / length(snp)) )
    # Find the mean
    mu.marker.miss <- mean(marker.miss, na.rm = T)

    ## Marker heterozygosity
    marker.het <- apply(X = geno.in, MARGIN = 2, FUN = function(snp) sum(snp == 0, na.rm = T) / sum(!is.na(snp)) )
    # Find the mean
    mu.marker.het <- mean(marker.het, na.rm = T)

    stat.list <- list( mu.MAF = mu.MAF, mu.marker.het = mu.marker.het, mu.marker.miss = mu.marker.miss,
                       MAF = MAF, marker.het = marker.het, marker.miss = marker.miss )

    return(stat.list)
  }

  # Entry stat calculator
  entry.stats <- function(geno.in) {
    # Calculate missing data proportion
    entry.miss <- apply(X = geno.in, MARGIN = 1, FUN = function(entry) return(sum(is.na(entry)) / length(entry)) )
    # Find the mean
    mu.entry.miss <- mean(entry.miss, na.rm = T)

    ## Entry heterozygosity
    entry.het <- apply(X = geno.in, MARGIN = 1, FUN = function(entry) sum(entry == 0, na.rm = T) / sum(!is.na(entry)) )
    # Find the mean
    mu.entry.het <- mean(entry.het, na.rm = T)

    # Output list
    stat.list <- list( mu.entry.het = mu.entry.het, mu.entry.miss = mu.entry.miss,
                       entry.het = entry.het, entry.miss = entry.miss )
    return(stat.list)
  }

  ## Marker stats
  ## Most marker stats calculate proportions by excluding missing data, except of course the
  # missing data calculator
  m.stats <- marker.stats(geno.in = geno.in)
  e.stats <- entry.stats(geno.in = geno.in)

  # Create a genotype table to output
  g.table <- table(geno.in)/(nrow(geno.in)*ncol(geno.in))

  ## Assemble lists for output
  output <- list( marker.stats = m.stats, entry.stats = e.stats, g.table = g.table)

  # Output graphs if told to
  if (graphs) {
    # Plot
    par(mfrow = c(2,3))
    # Site-frequency spectrum
    hist(m.stats$MAF, main = "Site Frequency Spectrum", xlab = "Minor Allele Frequency", xlim = c(0,0.5))
    # Missing data per marker
    hist(m.stats$marker.miss, main = "Frequency of Missing Marker Data", xlab = "Marker Missing Data Proportion", xlim = c(0,1))
    # Marker heterozygosity
    hist(m.stats$marker.het, main = "Frequency of Marker Heterozygosity", xlab = "Marker Heterozygosity", xlim = c(0,1))
    # Blank plot to take up space
    plot.new()
    # Missing data per line
    hist(e.stats$entry.miss, main = "Frequency of Missing Entry Data", xlab = "Line Missing Data Proportion", xlim = c(0,1))
    # Line heterozygosity
    hist(e.stats$entry.het, main = "Frequency of Entry Heterozygosity", xlab = "Entry Heterozygosity Proportion", xlim = c(0,1))
  }

  # Return the list
  return(output)

} # Close the function


### plot.PCA
##  This code takes a genetic relationship matrix and plots the principal components

plot.PCA <- function(Gmat = NULL, # A genetic relatioship matrix of dimensions n x n
                   x.PC = 1, # The number of the PC for the x-axis
                   y.PC = 2, # The number of the PC for the y-axis
                   color.factor = NULL, # A factor of categories by which to color the point. Must be of length n.
                   scree = FALSE, # Logical. If TRUE, display a scree plot
                   return.data = FALSE ) { # Logical. If true, the SVD results are returned

  # Error checking
  if (is.null(Gmat)) stop("No input given for 'Gmat'.")
  if (!is.matrix(Gmat)) { stop("The 'Gmat' input is not of class: matrix.") }
  if (nrow(Gmat) != ncol(Gmat)) stop("'Gmat' is not a square matrix.")
  if (is.null(color.factor)) {
    color.factor <- factor(rep(1, nrow(Gmat)))
  } else {
    if (!is.factor(color.factor)) color.factor <- factor(color.factor)
    if (length(color.factor) != nrow(Gmat)) stop("The length of 'color.factor' is not the same as the number of row/columns of 'Gmat'.")
  }

  # Perform the svd
  svd.out <- svd(Gmat)

  # Retrieve lambda
  lambda <- svd.out$d
  # Calculation the proportion of variation explained by each PC
  prop.var <- lambda / sum(lambda)

  # If scree is true, display the plot
  if (scree) {
    plot(x = 1:length(lambda),
         y = prop.var,
         xlab = "Principal Component",
         ylab = "Variation Explained by PC")

    # End the function
    return()
  }

  # Subset the principal components
  PC.x <- svd.out$u[,x.PC]
  PC.y <- svd.out$u[,y.PC]

  # Create pretty axis limits
  range.y <- range(pretty(range(PC.y)))
  ylim <- c( min(range.y), (max(range.y) + ( diff(range.y) / 4 )) )
  range.x <- range(pretty(range(PC.x)))
  xlim <- c( (min(range.x) - ( diff(range.x) / 6 )), max(range.x) )

  # Plot the PCA
  plot(x = PC.x,
       y = PC.y,
       xlab = paste("PC", x.PC, " (", round(prop.var[x.PC], 3) * 100, "%)", sep = ""),
       ylab = paste("PC", y.PC, " (", round(prop.var[y.PC], 3) * 100, "%)", sep = ""),
       ylim = ylim,
       xlim = xlim,
       pch = 16,
       col = color.factor)
  # Add a legend
  legend("topleft", legend = unique(color.factor), col = as.numeric(unique(color.factor)), pch = 16)

  # Output data if told
  if (return.data) return(svd.out)

} # Close the function
