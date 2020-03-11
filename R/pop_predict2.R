#' Predict genetic variance and genetic correlations using a deterministic model
#'
#' @description
#'
#' @param G.in See \code{G.in} in \code{\link[PopVar]{pop.predict}}.
#' @param y.in See \code{y.in} in \code{\link[PopVar]{pop.predict}}.
#' @param map.in See \code{map.in} in \code{\link[PopVar]{pop.predict}}.
#' @param crossing.table See \code{crossing.table} in \code{\link[PopVar]{pop.predict}}.
#' @param parents See \code{parents} in \code{\link[PopVar]{pop.predict}}.
#' @param tail.p See \code{tail.p} in \code{\link[PopVar]{pop.predict}}.
#' @param self.gen The number of selfing generations in the potential cross. Can be an integer or \code{Inf} for
#' recombinant inbreds. Note: \code{self.gen = 1} corresponds to an F2 population.
#' @param DH Indicator if doubled-haploids are to be induced after the number of selfing generations indicated by
#' \code{self.gen}. For example, if \code{self.gen = 0} and \code{DH = TRUE}, then doubled-haploids are asssumed
#' to be induced using gametes from F1 plants.
#' @param model See \code{models} in \code{\link[PopVar]{pop.predict}}. Only 1 model is allowed.
#'
#' @examples
#'
#' # Load data
#' data("phenos")
#' data("genos")
#' data("map")
#'
#' # Create 10, 2-way parent combinations
#' crosses <- as.data.frame(
#'    matrix(data = sample(row.names(genos), 20), nrow = 10, byrow = TRUE,
#'           dimnames = list(NULL, paste0("parent", 1:2))),
#'    stringsAsFactors = FALSE)
#'
#' # Format the genotype data
#' G_in <- as.data.frame(cbind( c("", row.names(genos)), rbind(colnames(genos), genos)) )
#'
#' # Run predictions
#' pred_out <- pop.predict2(G.in = G_in, y.in = phenos, map.in = map,
#'                          crossing.table = crosses)
#'
#'
#' @importFrom qtl mf.h
#' @importFrom rrBLUP mixed.solve
#'
#' @export
#'
pop.predict2 <- function(G.in, y.in, map.in, crossing.table, parents, tail.p = 0.1, self.gen = Inf, DH = FALSE,
                         model = c("rrBLUP", "BayesC")) {


  ## Check classes
  stopifnot(is.data.frame(G.in))
  stopifnot(is.data.frame(y.in))
  stopifnot(is.data.frame(map.in))

  # Check the number of markers in the map and the G.in objects
  if (ncol(G.in) - 1 != nrow(map.in))
    stop("The number of markers in G.in does not match the markers in map.in")

  ## Make sure there is no missing marker data
  if (any(is.na(G.in))) stop ("There must be no missing marker data.")

  # Self gen cannot be negative
  stopifnot(self.gen >= 0)

  # DH must be logical
  stopifnot(is.logical(DH))

  ## Make sure markers in G.in are in map.in
  markers_Gin <- as.character(unlist(G.in[1,-1, drop = TRUE]))
  markers_mapin <- map.in[,1, drop = TRUE]

  if (any(! markers_Gin %in% markers_mapin)) stop("The marker names in G.in are not all in map.in.")


  ## Filter map.in to remove markers with unknown position
  map.in[[1]] <- as.character(map.in[[1]])
  map.in[[3]] <- as.numeric(map.in[[3]])
  map.in_use <- subset(map.in, !is.na(map.in[[3]]))
  # Reorder based on chrom and then pos
  map.in_use <- map.in_use[order(map.in_use[[2]], map.in_use[[3]]), , drop = FALSE]

  # Get the names of the markers in the new map
  markers_mapin <- as.character(map.in_use[[1]])

  ## Subset G.in for markers in map.in_use
  G.in_use <- G.in[, c(1, which(markers_Gin %in% markers_mapin) + 1), drop = FALSE]



  # If the crossing table is not missing, check that the parents are in the G.in input
  if (!missing(crossing.table)) {
    parents <- unique(unlist(crossing.table))

    # Make sure the parent names are not factors
    crossing.table <- as.data.frame(sapply(X = crossing.table, as.character), stringsAsFactors = FALSE)

  } else {
    if (missing(parents))
      stop("If no crossing.table is provided, a list of parents must be supplied.")

    parents <- sort(parents)
    # Create a crossing table with all possible parent combinations
    crossing.table <- as.data.frame(t(combn(x = parents, m = 2)), stringsAsFactors = FALSE)
    names(crossing.table) <- c("parent1", "parent2")

  }

  if (any(!parents %in% G.in_use[,1,drop = T]))
    stop("Parents are not in G.in.")

  ## If self.gen is Inf and DH is T, error
  if (is.infinite(self.gen) & DH) stop("Infinite selfing generations and doubled-haploid production cannot both occur.")

  ## Set the factors of line names in y.in to those in the marker df and those in the y.in
  lines_G.in <- as.character(G.in_use[-1,1, drop = TRUE])
  # The levels of the y.in_use geno factor should be the entry names in G.in
  y.in[[1]] <- factor(x = y.in[[1]], levels = sort(lines_G.in))
  # Subset out NAs
  y.in_use <- subset(y.in, !is.na(y.in[[1]]))
  # Reorder lines
  y.in_use <- y.in_use[order(y.in_use[[1]]),, drop = FALSE]


  # Match arguments
  model <- match.arg(model)

  ## Number of traits and trait names
  n_traits <- ncol(y.in_use) - 1
  trait_names <- colnames(y.in_use)[-1]




  ## Calculate the expected genetic variance and covariance of markers

  ## Create an empty matrix
  marker_names <- markers_mapin
  covar <- matrix(0, nrow = nrow(map.in_use), ncol = nrow(map.in_use), dimnames = list(marker_names, marker_names))

  # Split markers by chromosome
  map.in.chr <- split(map.in_use, map.in_use[,2, drop = FALSE])
  markers_chr <- lapply(map.in.chr, "[[", 1)

  # Calculate separate centimorgan distance matrices per chromosome
  chr_cM <- lapply(X = map.in.chr, FUN = function(x) as.matrix(dist(x[,3,drop = FALSE])))
  # Convert to recombination distance
  chr_c <- lapply(X = chr_cM, FUN = mf.h)



  # Calculate the LD covariance

  ## This depends on the number of selfing generations
  ## Formulae are taken from Leheremeir et al 2017

  # If DH and self.gen = 0, DH's are formed from the F1
  if (DH & self.gen == 0) {
    chr_covar <- lapply(X = chr_c, FUN = function(c) 1 - (2 * c))

  } else if (DH & self.gen > 0) {
    covar_selfing <- lapply(X = chr_c, FUN = function(c) {
      base <- 0.5 * (1 - (2 * c))
      Reduce(f = `+`, x = lapply(X = seq(self.gen), FUN = function(k) base ^ k)) })

    covar_final <- lapply(X = chr_c, FUN = function(c) (0.5 * (1 - (2 * c))) ^ self.gen)

    # Sum
    chr_covar <- mapply(FUN = `+`, covar_selfing, covar_final)

  } else if (!DH & is.finite(self.gen)) {
    chr_covar <- lapply(X = chr_c, FUN = function(c) {
      base <- 0.5 * (1 - (2 * c))
      Reduce(f = `+`, x = lapply(X = seq(self.gen), FUN = function(k) base ^ k)) })

  } else if (!DH & is.infinite(self.gen)) {
    chr_covar <- lapply(X = chr_c, FUN = function(c) (1 - (2 * c)) / (1 + (2 * c)))

  }



  ## Convert the marker matrix into a useable matrix form
  G.in_pred <- sapply(X = G.in_use[-1, -1], function(x) as.numeric(as.character(x)))
  row.names(G.in_pred) <- as.character(G.in_use[-1,1])
  colnames(G.in_pred) <- as.character(unlist(G.in_use[1,-1]))
  # Reorder markers and lines
  G.in_pred <- G.in_pred[order(row.names(G.in_pred)), marker_names]


  ## Create a model.frame from the phenotypes; extract the column name of genotypes/lines
  geno_colname <- colnames(y.in_use)[1]

  # Subset using factors in y.in
  # Model matrix for genotypes
  Zg <- model.matrix(object = as.formula(paste0("~ -1 + ", geno_colname)), data = y.in_use)
  row.names(Zg) <- y.in_use[[1]]
  M <- Zg %*% G.in_pred


  # Make sure the markers are coded correlated
  if (!all(M %in% c(-1, 1))) stop("All markers must be coded as -1 or 1. Heterzygotes are not allowed.")


  ## Calculate marker effects for each trait
  if (model == "rrBLUP") {

    ## Apply a function over each trait
    marker_effect_out <- lapply(X = y.in_use[-1], FUN = function(y) {

      ## Solve the mixed model
      fit <- mixed.solve(y = y, Z = M, method = "REML")

      # Return marker effects and the grand mean
      list(effects = as.matrix(fit$u), grand_mean = fit$beta)

    })


    ## Create a complete matrix of marker effects for each trait
    mar_eff_mat <- do.call("cbind", lapply(marker_effect_out, "[[", "effects"))
    mar_eff_mat <- structure(mar_eff_mat, dimnames = list(row.names(mar_eff_mat), names(marker_effect_out)))

    mar_beta_mat <- do.call("cbind", lapply(marker_effect_out, "[[", "grand_mean"))
    mar_beta_mat <- structure(mar_beta_mat, dimnames = list(row.names(mar_beta_mat), names(marker_effect_out)))


  } else if (model == "BayesC") {
    stop("Other models not supported.")

  } else {
    stop("Other models not supported.")
  }


  ## Predicted genotypic value of all genotypes + grand mean
  pgvs <- (M %*% mar_eff_mat) + matrix(mar_beta_mat, ncol = ncol(mar_eff_mat), nrow = nrow(M), byrow = TRUE)
  row.names(pgvs) <- y.in_use[[1]]


  ## Calculate the pairwise product of marker effects for each trait, separated by chromosome
  ## Then multiply by the LD covariance
  intra_trait_covar <- list()

  # Iterate over traits
  for (i in seq(ncol(mar_eff_mat))) {

    # Split the marker effect matrix by chromosome and calculate pairwise product
    mar_eff_prod_chr <- lapply(X = markers_chr, function(marks) tcrossprod(mar_eff_mat[marks, i, drop = F]))

    # The covariance is the QTL effect product multiplied by the expected D
    intra_trait_covar[[trait_names[i]]] <- mapply(mar_eff_prod_chr, chr_covar, FUN = `*`, SIMPLIFY = FALSE)

  }


  if (n_traits > 1) {

    ## Calculate the pairwise product of marker effects for each pair of traits
    # First create a vector of trait index combinations
    trait_ind_combn <- t(combn(x = seq(n_traits), m = 2, simplify = TRUE))
    trait_combn_name <- combn(x = trait_names, m = 2, paste0, collapse = ":")

    # Empty matrix for correlations
    trait_corG_mat <- matrix(data = NA, nrow = n_traits, ncol = n_traits, dimnames = list(trait_names, paste0("cor_w_", trait_names)))
    # Convert to a distance version for easy addition
    trait_corG_dist <- as.dist(trait_corG_mat)

    # List of inter-trait covariances
    inter_trait_covar <- list()

    # Iterate over trait indices combinations
    for (i in  seq(nrow(trait_ind_combn))) {

      # Get the marker effects for that combination
      mar_eff_combn <- mar_eff_mat[,trait_ind_combn[i,], drop = FALSE]

      # Calculate the pairwise product by chromosome
      mar_eff_prod_chr <- lapply(X = markers_chr, function(marks) tcrossprod(mar_eff_combn[marks,1,drop = FALSE], mar_eff_combn[marks,2,drop = FALSE]))

      # The covariance is the QTL effect product multiplied by the expected D
      inter_trait_covar[[trait_combn_name[i]]] <- mapply(mar_eff_prod_chr, chr_covar, FUN = `*`, SIMPLIFY = FALSE)

    }

  } else {
    inter_trait_covar <- NULL

  }


  # Determine the k_sp from the the tail.p
  k_sp <- mean(qnorm(p = seq((1 - tail.p), 0.999999, 0.000001)))


  ## Create a list to store dfs
  cross_predictions <- vector("list", length = nrow(crossing.table))

  ## Iterate over the pairs of parents
  for (j in seq(nrow(crossing.table))) {

    # Character vector of the two parents
    pars <- as.character(crossing.table[j,1:2])

    ## Subset the genotype matrix using the parents
    par_geno <- G.in_pred[pars,,drop = FALSE]

    # Which markers are segregating?
    mar_seg <- names(which(colMeans(par_geno) == 0))
    # Split by chromsome
    mar_seg_chr <- lapply(markers_chr, intersect, mar_seg)

    # Find the parent 1 genotype of those markers and take the crossproduct for multiplication
    par1_mar_seg <- crossprod(par_geno[1,mar_seg, drop = FALSE])
    # Split by chromosome
    par1_mar_seg_chr <- lapply(mar_seg_chr, function(marks) par1_mar_seg[marks, marks])


    ## Predictions
    # Cross mean
    pred_mu_j <- colMeans(pgvs[pars,,drop = FALSE])

    # Genetic variance - use the pred_mu_j vector as a template
    pred_varG_j <- pred_mu_j

    # Iterate over traits
    # Note that the QTL covariance matrix includes the variance of each QTL on the diagonal, so the sum of the matrix
    # is the variance + 2 * covariance
    for (i in seq(n_traits)) {
      pred_varG_j[i] <- sum(mapply(par1_mar_seg_chr, intra_trait_covar[[i]], FUN = function(x, y) sum(x * y[colnames(x), colnames(x)])))
    }


    # Genetic correlations between traits, if more than one trait
    if (!is.null(inter_trait_covar)) {

      ## Iterate over trait combinations
      for (i in seq_along(inter_trait_covar)) {

        ## Calculate the covariance between the pair of traits
        trait_pair_cov <- sum(mapply(par1_mar_seg_chr, inter_trait_covar[[i]], FUN = function(x, y) sum(x * y[colnames(x), colnames(x)])))

        # Subset predicted variance
        trait_pair_varG <- pred_varG_j[trait_ind_combn[i,]]

        # Calculate correlation and save
        trait_corG_dist[i] <- trait_pair_cov / prod(sqrt(trait_pair_varG))

      }

      ## Convert distance object to matrix
      pred_corG_mat <- as.matrix(trait_corG_dist)
      dimnames(pred_corG_mat) <- dimnames(trait_corG_mat)
      diag(pred_corG_mat) <- NA

      ## Calculate correlated progeny mean
      response_trait_varG <- matrix(pred_varG_j, nrow = length(pred_varG_j), ncol = length(pred_varG_j), byrow = TRUE)
      correlated_response <- k_sp * pred_corG_mat * sqrt(response_trait_varG)
      pred_mu_j_mat <- matrix(pred_mu_j, nrow = length(pred_mu_j), ncol = length(pred_mu_j), byrow = TRUE)
      pred_cor_musp_low <- pred_mu_j_mat - correlated_response
      pred_cor_musp_high <- pred_mu_j_mat + correlated_response

      # Change names
      colnames(pred_cor_musp_low) <- paste0("pred_cor_musp_low_", trait_names)
      colnames(pred_cor_musp_high) <- paste0("pred_cor_musp_high_", trait_names)

    } else {
      pred_corG_mat <- pred_cor_musp_low <- pred_cor_musp_high <- NULL

    }


    ## Save the results as a data.frame
    cross_predictions[[j]] <- data.frame(parent1 = pars[1], parent2 = pars[2], trait = trait_names,
                                         cbind(pred_mu = pred_mu_j, pred_varG = pred_varG_j, pred_corG_mat, pred_cor_musp_low, pred_cor_musp_high),
                                         stringsAsFactors = FALSE, row.names = NULL)

  } # End loop

  ## Bind rows
  cross_predictions1 <- do.call("rbind", cross_predictions)

  ## Calculate response predictions (superior progeny, correlated response, etc.)
  pred_response <- (k_sp * sqrt(cross_predictions1$pred_varG))

  # Superior progeny mean
  cross_predictions1[["pred_musp_low"]] <- cross_predictions1$pred_mu - pred_response
  cross_predictions1[["pred_musp_high"]] <- cross_predictions1$pred_mu + pred_response

  ## Re-order columns
  cor_W_cols <- grep(pattern = paste0(trait_names, collapse = "|"), x = names(cross_predictions1))
  cross_predictions2 <- cross_predictions1[, c(setdiff(seq(ncol(cross_predictions1)), cor_W_cols), cor_W_cols)]

  ## Return the predictions
  return(cross_predictions2)

} # Close function




#' Predict genetic variance and genetic correlations using a deterministic model
#'
#' @describeIn pop.predict2
#'
#'
#' @param M A Matrix of marker genotypes of dimensions \code{nLine} x \code{nMarker}, coded as
#' -1, 0, and 1.
#' @param y.in A data frame of phenotypic means. The first column should include the entry name and
#' subsequent columns should include phenotypic values. Ignored if \code{marker.effects} is passed.
#' @param marker.effects A data frame of marker effects. The first column should include the marker name and
#' subsequent columns should include the marker effects.
#' @param map.in See \code{map.in} in \code{\link[PopVar]{pop.predict}}.
#' @param crossing.table See \code{crossing.table} in \code{\link[PopVar]{pop.predict}}.
#' @param tail.p See \code{tail.p} in \code{\link[PopVar]{pop.predict}}.
#' @param self.gen The number of selfing generations in the potential cross. Can be an integer or \code{Inf} for
#' recombinant inbreds. Note: \code{self.gen = 1} corresponds to an F2 population.
#' @param DH Indicator if doubled-haploids are to be induced after the number of selfing generations indicated by
#' \code{self.gen}. For example, if \code{self.gen = 0} and \code{DH = TRUE}, then doubled-haploids are asssumed
#' to be induced using gametes from F1 plants.
#' @param model See \code{models} in \code{\link[PopVar]{pop.predict}}. Only 1 model is allowed.
#' @param n.core Number of cores for parallelization; only works on a Linux or Mac OS operating system.
#'
#' @examples
#'
#' # Load data
#' data("phenos")
#' data("genos")
#' data("map")
#'
#' # Create 10, 2-way parent combinations
#' crosses <- as.data.frame(
#'    matrix(data = sample(row.names(genos), 20), nrow = 10, byrow = TRUE,
#'           dimnames = list(NULL, paste0("parent", 1:2))),
#'    stringsAsFactors = FALSE)
#'
#' # Run predictions
#' pred_out <- pop_predict2(M = genos, y.in = phenos, map.in = map,
#'                          crossing.table = crosses)
#'
#'
#' ## Pass marker effects instead of phenotypes
#' # First calculate marker effects
#' phenos2 <- as.matrix(phenos[,-1]); row.names(phenos2) <- phenos[,1]
#' phenos2 <- phenos2[row.names(genos),]
#'
#' mar_eff <- apply(X = phenos2, MARGIN = 2, FUN = function(y) mixed.solve(y = y, Z = genos)$u)
#' marker_effects <- data.frame(marker = row.names(mar_eff), mar_eff, stringsAsFactors = FALSE)
#'
#' pred_out <- pop_predict2(M = genos, marker.effects = marker_effects, map.in = map,
#'                          crossing.table = crosses, self.gen = 6)
#'
#'
#'
#' @importFrom qtl mf.h
#' @importFrom rrBLUP mixed.solve
#' @importFrom parallel mclapply
#'
#' @export
#'
pop_predict2 <- function(M, y.in, marker.effects, map.in, crossing.table, tail.p = 0.1,
                         self.gen = Inf, DH = FALSE, model = c("rrBLUP", "BayesC"), n.core = 1) {


  ###################
  # Error handling

  ## Check classes
  stopifnot(is.matrix(M))
  stopifnot(is.data.frame(map.in))

  # Check the number of markers in the map and the G.in objects
  if (ncol(M) != nrow(map.in))
    stop("The number of markers in G.in does not match the markers in map.in")

  ## Make sure there is no missing marker data
  if (any(is.na(M))) stop ("There must be no missing marker data.")

  # Self gen cannot be negative
  stopifnot(self.gen >= 0)
  # DH must be logical
  stopifnot(is.logical(DH))

  ## Make sure markers in G.in are in map.in
  stopifnot(!is.null(colnames(M)))
  stopifnot(!is.null(row.names(M)))

  # Make sure the markers are coded correctly
  if (!all(M %in% c(-1, 1))) stop("All markers must be coded as -1 or 1. Heterzygotes are not allowed.")


  markers_M <- colnames(M)
  markers_mapin <- map.in[[1]]

  if (any(! markers_M %in% markers_mapin)) stop("The marker names in M are not all in map.in.")

  # Get the names of genotyped entries
  geno_lines <- row.names(M)


  ## Make sure one of y.in or marker.effects are not missing
  if (missing(y.in) & missing(marker.effects))
    stop("You must pass one of 'y.in' or 'marker.effects.'")

  ## Error check depending on what is passed
  if (!missing(marker.effects)) {

    # Check markers
    markers_maref <- marker.effects[[1]]

    if (! all(markers_M %in% markers_maref) ) stop("The marker names in M are not all in marker.effects")

    ## Number of traits and trait names
    n_traits <- ncol(marker.effects) - 1
    trait_names <- colnames(marker.effects)[-1]

    # Set boolean for later
    calc_marker_eff <- FALSE

    # Else check y.in
  } else {

    stopifnot(is.data.frame(y.in))

    # All phenotyped lines should be genotyped
    if (! all(y.in[[1]] %in% geno_lines) ) stop("All entries in 'y.in' should have marker data in 'M'.")

    ## Set the factors of line names in y.in to those in the marker df and those in the y.in
    y.in[[1]] <- factor(x = y.in[[1]], levels = sort(geno_lines))
    # Subset out NAs
    y.in_use <- subset(y.in, !is.na(y.in[[1]]))
    # Reorder lines
    y.in_use <- y.in_use[order(y.in_use[[1]]),, drop = FALSE]


    ## Number of traits and trait names
    n_traits <- ncol(y.in_use) - 1
    trait_names <- colnames(y.in_use)[-1]

    # Set boolean for later
    calc_marker_eff <- TRUE

  }

  # Match arguments
  model <- match.arg(model)



  # Reorder map based on chrom and then pos
  map.in_use <- map.in[order(map.in[[2]], map.in[[3]]), , drop = FALSE]

  # Get the names of the markers in the new map
  markers_mapin <- as.character(map.in_use[[1]])


  # Make sure the parent names are not factors
  crossing.table <- as.data.frame(sapply(X = crossing.table, as.character), stringsAsFactors = FALSE)

  ## Get list of unique parents from the crossing.table
  parents <- unique(unlist(crossing.table))


  if (!all(parents %in% row.names(M)))
    stop("Parents are not in G.in.")


  ## Fit models to calculate marker effects, if necessary
  if (calc_marker_eff) {

    ## Create a model.frame from the phenotypes; extract the column name of genotypes/lines
    geno_colname <- colnames(y.in_use)[1]

    # Subset using factors in y.in
    # Model matrix for genotypes
    Zg <- model.matrix(object = as.formula(paste0("~ -1 + ", geno_colname)), data = y.in_use)
    row.names(Zg) <- y.in_use[[1]]
    M1 <- Zg %*% M[, markers_mapin, drop = FALSE]

    ## Calculate marker effects for each trait
    if (model == "rrBLUP") {

      ## Apply a function over each trait
      marker_effect_out <- lapply(X = y.in_use[-1], FUN = function(y) {

        ## Solve the mixed model
        fit <- mixed.solve(y = y, Z = M1, method = "REML")

        # Return marker effects and the grand mean
        list(effects = as.matrix(fit$u), grand_mean = fit$beta)

      })


      ## Create a complete matrix of marker effects for each trait
      mar_eff_mat <- do.call("cbind", lapply(marker_effect_out, "[[", "effects"))
      mar_eff_mat <- structure(mar_eff_mat, dimnames = list(row.names(mar_eff_mat), names(marker_effect_out)))

      mar_beta_mat <- do.call("cbind", lapply(marker_effect_out, "[[", "grand_mean"))
      mar_beta_mat <- structure(mar_beta_mat, dimnames = list(row.names(mar_beta_mat), names(marker_effect_out)))


    } else if (model == "BayesC") {
      stop("Other models not supported.")

    } else {
      stop("Other models not supported.")
    }

    # Else create a matrix of marker ordered marker effects
  } else {

    # Create matrix
    mar_eff_mat <- as.matrix(marker.effects[,-1,drop = FALSE])
    row.names(mar_eff_mat) <- marker.effects[[1]]
    mar_eff_mat <- mar_eff_mat[markers_mapin,,drop = FALSE]

    # Set the grand mean to zero
    mar_beta_mat <- matrix(0, nrow = 1, ncol = ncol(mar_eff_mat),
                           dimnames = list(NULL, colnames(mar_eff_mat)))

  }

  # Reorder markers
  M1 <- M[, markers_mapin, drop = FALSE]

  ## Create an empty matrix
  marker_names <- markers_mapin
  covar <- matrix(0, nrow = nrow(map.in_use), ncol = nrow(map.in_use), dimnames = list(marker_names, marker_names))

  # Split markers by chromosome
  map.in.chr <- split(map.in_use, map.in_use[,2, drop = FALSE])
  markers_chr <- lapply(map.in.chr, "[[", 1)

  # Calculate separate centimorgan distance matrices per chromosome
  chr_cM <- lapply(X = map.in.chr, FUN = function(x) as.matrix(dist(x[,3,drop = FALSE])))
  # Convert to recombination distance
  chr_c <- lapply(X = chr_cM, FUN = mf.h)



  # Calculate the LD covariance

  ## This depends on the number of selfing generations
  ## Formulae are taken from Leheremeir et al 2017

  # If DH and self.gen = 0, DH's are formed from the F1
  if (DH & self.gen == 0) {
    chr_covar <- lapply(X = chr_c, FUN = function(c) 1 - (2 * c))

  } else if (DH & self.gen > 0) {
    covar_selfing <- lapply(X = chr_c, FUN = function(c) {
      base <- 0.5 * (1 - (2 * c))
      Reduce(f = `+`, x = lapply(X = seq(self.gen), FUN = function(k) base ^ k)) })

    covar_final <- lapply(X = chr_c, FUN = function(c) (0.5 * (1 - (2 * c))) ^ self.gen)

    # Sum
    chr_covar <- mapply(FUN = `+`, covar_selfing, covar_final)

  } else if (!DH & is.finite(self.gen)) {
    chr_covar <- lapply(X = chr_c, FUN = function(c) {
      base <- 0.5 * (1 - (2 * c))
      Reduce(f = `+`, x = lapply(X = seq(self.gen), FUN = function(k) base ^ k)) })

  } else if (!DH & is.infinite(self.gen)) {
    chr_covar <- lapply(X = chr_c, FUN = function(c) (1 - (2 * c)) / (1 + (2 * c)))

  }


  ## Predicted genotypic value of all genotypes + grand mean
  pgvs <- (M %*% mar_eff_mat) + matrix(mar_beta_mat, ncol = ncol(mar_eff_mat), nrow = nrow(M), byrow = TRUE)


  ## Calculate the pairwise product of marker effects for each trait, separated by chromosome
  ## Then multiply by the LD covariance
  intra_trait_covar <- list()

  # Iterate over traits
  for (i in seq(ncol(mar_eff_mat))) {

    # Split the marker effect matrix by chromosome and calculate pairwise product
    mar_eff_prod_chr <- lapply(X = markers_chr, function(marks) tcrossprod(mar_eff_mat[marks, i, drop = F]))

    # The covariance is the QTL effect product multiplied by the expected D
    intra_trait_covar[[trait_names[i]]] <- mapply(mar_eff_prod_chr, chr_covar, FUN = `*`, SIMPLIFY = FALSE)

  }


  if (n_traits > 1) {

    ## Calculate the pairwise product of marker effects for each pair of traits
    # First create a vector of trait index combinations
    trait_ind_combn <- t(combn(x = seq(n_traits), m = 2, simplify = TRUE))
    trait_combn_name <- combn(x = trait_names, m = 2, paste0, collapse = ":")

    # Empty matrix for correlations
    trait_corG_mat <- matrix(data = NA, nrow = n_traits, ncol = n_traits, dimnames = list(trait_names, paste0("cor_w_", trait_names)))
    # Convert to a distance version for easy addition
    trait_corG_dist <- as.dist(trait_corG_mat)

    # List of inter-trait covariances
    inter_trait_covar <- list()

    # Iterate over trait indices combinations
    for (i in  seq(nrow(trait_ind_combn))) {

      # Get the marker effects for that combination
      mar_eff_combn <- mar_eff_mat[,trait_ind_combn[i,], drop = FALSE]

      # Calculate the pairwise product by chromosome
      mar_eff_prod_chr <- lapply(X = markers_chr, function(marks) tcrossprod(mar_eff_combn[marks,1,drop = FALSE], mar_eff_combn[marks,2,drop = FALSE]))

      # The covariance is the QTL effect product multiplied by the expected D
      inter_trait_covar[[trait_combn_name[i]]] <- mapply(mar_eff_prod_chr, chr_covar, FUN = `*`, SIMPLIFY = FALSE)

    }

  } else {
    inter_trait_covar <- NULL

  }


  # Determine the k_sp from the the tail.p
  k_sp <- mean(qnorm(p = seq((1 - tail.p), 0.999999, 0.000001)))


  ## Create a list to store dfs
  cross_predictions <- vector("list", length = nrow(crossing.table))

  ## Iterate over the pairs of parents
  for (j in seq(nrow(crossing.table))) {

    # Character vector of the two parents
    pars <- as.character(crossing.table[j,1:2])

    ## Subset the genotype matrix using the parents
    par_geno <- M[pars,,drop = FALSE]

    # Which markers are segregating?
    mar_seg <- names(which(colMeans(par_geno) == 0))
    # Split by chromsome
    mar_seg_chr <- lapply(markers_chr, intersect, mar_seg)

    # Find the parent 1 genotype of those markers and take the crossproduct for multiplication
    par1_mar_seg <- crossprod(par_geno[1,mar_seg, drop = FALSE])
    # Split by chromosome
    par1_mar_seg_chr <- lapply(mar_seg_chr, function(marks) par1_mar_seg[marks, marks])


    ## Predictions
    # Cross mean
    pred_mu_j <- colMeans(pgvs[pars,,drop = FALSE])

    # Genetic variance - use the pred_mu_j vector as a template
    pred_varG_j <- pred_mu_j

    # Iterate over traits
    # Note that the QTL covariance matrix includes the variance of each QTL on the diagonal, so the sum of the matrix
    # is the variance + 2 * covariance
    for (i in seq(n_traits)) {
      pred_varG_j[i] <- sum(mapply(par1_mar_seg_chr, intra_trait_covar[[i]], FUN = function(x, y) sum(x * y[colnames(x), colnames(x)])))
    }


    # Genetic correlations between traits, if more than one trait
    if (!is.null(inter_trait_covar)) {

      ## Iterate over trait combinations
      for (i in seq_along(inter_trait_covar)) {

        ## Calculate the covariance between the pair of traits
        trait_pair_cov <- sum(mapply(par1_mar_seg_chr, inter_trait_covar[[i]], FUN = function(x, y) sum(x * y[colnames(x), colnames(x)])))

        # Subset predicted variance
        trait_pair_varG <- pred_varG_j[trait_ind_combn[i,]]

        # Calculate correlation and save
        trait_corG_dist[i] <- trait_pair_cov / prod(sqrt(trait_pair_varG))

      }

      ## Convert distance object to matrix
      pred_corG_mat <- as.matrix(trait_corG_dist)
      dimnames(pred_corG_mat) <- dimnames(trait_corG_mat)
      diag(pred_corG_mat) <- NA

      ## Calculate correlated progeny mean
      response_trait_varG <- matrix(pred_varG_j, nrow = length(pred_varG_j), ncol = length(pred_varG_j), byrow = TRUE)
      correlated_response <- k_sp * pred_corG_mat * sqrt(response_trait_varG)
      pred_mu_j_mat <- matrix(pred_mu_j, nrow = length(pred_mu_j), ncol = length(pred_mu_j), byrow = TRUE)
      pred_cor_musp_low <- pred_mu_j_mat - correlated_response
      pred_cor_musp_high <- pred_mu_j_mat + correlated_response

      # Change names
      colnames(pred_cor_musp_low) <- paste0("pred_cor_musp_low_", trait_names)
      colnames(pred_cor_musp_high) <- paste0("pred_cor_musp_high_", trait_names)

    } else {
      pred_corG_mat <- pred_cor_musp_low <- pred_cor_musp_high <- NULL

    }


    ## Save the results as a data.frame
    cross_predictions[[j]] <- data.frame(parent1 = pars[1], parent2 = pars[2], trait = trait_names,
                                         cbind(pred_mu = pred_mu_j, pred_varG = pred_varG_j, pred_corG_mat, pred_cor_musp_low, pred_cor_musp_high),
                                         stringsAsFactors = FALSE, row.names = NULL)

  } # End loop

  ## Bind rows
  cross_predictions1 <- do.call("rbind", cross_predictions)

  ## Calculate response predictions (superior progeny, correlated response, etc.)
  pred_response <- (k_sp * sqrt(cross_predictions1$pred_varG))

  # Superior progeny mean
  cross_predictions1[["pred_musp_low"]] <- cross_predictions1$pred_mu - pred_response
  cross_predictions1[["pred_musp_high"]] <- cross_predictions1$pred_mu + pred_response

  ## Re-order columns
  cor_W_cols <- grep(pattern = paste0(trait_names, collapse = "|"), x = names(cross_predictions1))
  cross_predictions2 <- cross_predictions1[, c(setdiff(seq(ncol(cross_predictions1)), cor_W_cols), cor_W_cols)]

  ## Return the predictions
  return(cross_predictions2)

} # Close function
