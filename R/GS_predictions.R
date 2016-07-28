### G_predict.R function to predict the breeding value of a population using genome-wide markers
## Author: Jeff Neyhart
## University of Minnesota, St. Paul
## January 15, 2016

g.predict <- function(y.train = NULL, # A n x p matrix of training phenotype data, where rows are lines (with row names being the line names) and columns are phenotypic data for a trait or trait-environment combination (and column names refer to trait names or trait-environment combinations)
                      g.train = NULL, # A n x m matrix of training marker data, coded as z { -1, 0 , 1 } with row names being the line names and column names being the marker names
                      g.pred = NULL, # A n x m matrix of prediction marker data, coded as z { -1, 0 , 1 } with row names being the line names and column names being the marker names
                      Amat = NULL, # A n x n additive relationship matricies of training and prediction lines, with row names being line names
                      method = NULL, # A character vector of prediction methods to use. Options are "rrBLUP", "GBLUP", "BL", and "RKHS". If NULL, 10-fold cross validation is used to determine the best method per trait or trait-environment combination
                      n.iter = 1200, # Number of iterations for bayesian regression
                      n.burn = 200, # Number of iteration to burn in for bayesian regression
                      output = "csv", # Output format. Options are "csv" for a .csv file, or "rdata" for an .RData file
                      out.filename = NULL # The output basename, excluding the extension. If NULL, no file will be outputted
                      ) {
  
  # Load libraries
  library(rrBLUP)
  library(BGLR)
  
  # Available methods
  pred.methods <- c("rrBLUP", "GBLUP", "BL", "RKHS")

  # Error reporting
  if(is.null(y.train)) { stop("The y.train input was not specified") }
  if(all(is.null(g.train), is.null(g.pred), is.null(Amat))) { stop("Need to specify genotype data as either 1) marker matrices or 2) an additive relationship matrix") }
  if(!is.null(method) && !all(method %in% pred.methods)) { stop("One or more of the methods in the methods input is incorrect.") }
  if(!any(output == c("csv", "rdata"))) { stop("The output argument is not an option.") }
  
  # Convert to matricies
  y.train <- as.matrix(y.train)
  g.train <- as.matrix(g.train)
  g.pred <- as.matrix(g.pred)
  
  # Retrieve line names
  t.names <- list(row.names(y.train), row.names(g.train))
  # Find the shorter list and make that the final line list
  t.names <- t.names[[which.min(sapply(X = t.names, FUN = length))]]
  p.names <- row.names(g.pred)
  
  # Retrieve trait-env names
  te <- colnames(y.train)
  # Separate into traits and environments
  traits <- unique(sapply(X = strsplit(x = te, split = "_"), FUN = function(x) return(x[1])))
  envs <- unique(sapply(X = strsplit(x = te, split = "_"), FUN = function(x) return(paste(x[2], x[3], sep = ""))))
  envs <- gsub(pattern = "NA", replacement = "", x = envs)
  envs <- envs[envs != ""]
  
  # Subset the data based on line names
  y.train <- y.train[t.names,]
  g.train <- g.train[t.names,]
  
  # Sort the data
  y.train <- y.train[order(row.names(y.train)),]
  g.train <- g.train[order(row.names(y.train)),]
  g.pred <- g.pred[order(row.names(g.pred)),]
  
  # Create a relationship matrix
  g.tot <- rbind(g.train, g.pred)
  Amat.tot <- A.mat(X = g.tot, min.MAF = 0, max.missing = 1)
  
  # Count lines and markers
  n.t <- length(t.names)
  n.p <- length(p.names)
  m <- ncol(g.train)
  
  # Define functions
  # RR-BLUP marker effects
  rrblup.solve <- function(y, Z, K, X) {
    solve.out <- mixed.solve(y = y, Z = Z, K = K, X = X, SE = F, return.Hinv = F)
    out <- list(m.effects = as.matrix(solve.out$u), mu = solve.out$beta)
    return(out)
  }
  # GBLUP predictions
  gblup.solve <- function(y, Z, K, X) {
    solve.out <- mixed.solve(y = y, Z = Z, K = K, X = X, SE = F, return.Hinv = F)
    yhat <- solve.out$u + as.numeric(solve.out$beta)
    return(yhat)
  }
  # BL marker effects
  bl.solve <- function(y, Z, K, X) {
    ETA <- list(list(X = X, model = "FIXED"), list(X = Z, model = "BL", list(K = K, model = "RKHS")))
    solve.out <- BGLR(y = y, ETA = ETA, verbose = F, nIter = n.iter, burnIn = n.burn)  
    out <- list(m.effects = as.matrix(solve.out$ETA[[2]]$b, mu = solve.out$mu))
    return(out)
  }
  # RKHS marker effects
  rkhs.solve <- function(y, Z, K, X) {
    ETA <- list(list(X = X, model = "FIXED"), list(X = Z, model = "BL", list(K = K, model = "RKHS")))
    solve.out <- BGLR(y = y, ETA = ETA, verbose = F, nIter = n.iter, burnIn = n.burn)  
    return(solve.out$yhat)
  }
    
  # List of functions
  pred.funcs <- list(rrBLUP = rrblup.solve, GBLUP = gblup.solve, BL = bl.solve, RKHS = rkhs.solve)
  
  # Create a function to make predictions (for each trait) of a trait in each environment
  
  # Create output lists / data.frames
  pred.list <- list()
  
  # Assign method and design matricies
  pred.method = "rrBLUP"
  
  X = matrix(1, n.t)
  K = diag(m)
  
  # Make predictions
  # For each trait, create a matrix of values for just that trait
  ## Then for each trait in each environment, make predictions
  for (trait in traits) {
    # Notify
    cat(paste("Making predictions for ", trait, ".\n", sep = ""))
    # Create a results list
    trait.list <- list()
    # Find the index of the y.train data for that trait
    trait.index <- grepl(pattern = trait, x = colnames(y.train))
    # Subset the y.train matrix
    y.train.trait <- subset.matrix(x = y.train, select = trait.index)
    
    # For each column (env) in the subset, make predictions
    for (i in 1:ncol(y.train.trait)) {
      y <- as.numeric(y.train.trait[,i])
      env <- colnames(y.train.trait)[i]
      solve.out <- pred.funcs[[pred.method]](y = y, Z = g.train, K = K, X = X)
      GEBVs <- g.pred %*% solve.out$m.effects
      yhat <- GEBVs + as.numeric(solve.out$mu)
      # Add to the results list
      trait.list[[env]] <- yhat
    }
    
    # Make the list into a matrix
    trait.mat <- do.call("cbind", trait.list); colnames(trait.mat) <- names(trait.list)
    # Add to the pred.list
    pred.list[[trait]] <- trait.mat
    
  } # Close the per-trait loop
  
  # Make the list into a matrix
  pred.mat <- do.call("cbind", pred.list)
  
  # Output
  # Only output if told
  if(!is.null(out.filename)) {
    # Remove extension
    out.filename <- sub(pattern = "^([^.]*).*", replacement = "\\1", x = out.filename)
    if (output == "csv") {
      ext <- ".csv"
      filename <- paste(out.filename, ext, sep = "")
      write.csv(x = pred.mat, file = filename, quote = F)
    }
    if (output == "rdata") {
      ext <- ".RData"
      filename <- paste(out.filename, ext, sep = "")
      save(list = "pred.mat", file = filename)
    }
  }
  
  # Return the matrix
  return(pred.mat)
  
} # Close the function