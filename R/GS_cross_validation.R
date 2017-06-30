## Script to define some functions for genomic selection cross-validation

# Dependencies
# Load libraries
library(rrBLUP)
# library(BGLR)
# library(randomForest)

## Testing. DO NOT USE ####
# # Define the function for fractional cross-validation (i.e. 60:40, 50:50, etc)
# frac.xval <- function(g.in = NULL, # A n x m genotype matrix (coded as -1, 0, 1), with the row.names of g.in equal to the line names
#                       y.in = NULL, # A n x t matrix of phentypes, with the row.names of y.in equal to the line names
#                       frac.train = 0.70, # The fraction of the lines to use as the training set
#                       reps = 250, # The number of iterations of fractional cross-validation
#                       CV.nIter = 1500, # Number of interations of Gibbs sampling to perform using the Bayesian models
#                       CV.burnIn = 500, # Number of CV.nIter to discard before using the remainder to establish the prior distribution for Bayesian methods
#                       models = c("rrBLUP", "BayesA", "BayesB", "BayesC", "BL", "BRR", "RKHS"), # Vector of model names to employ
#                       n.cores = NULL # Optional number of cores to use when parallelizing
#                     ) {
#   
#   # Vector of all models available
#   all.models <- c("rrBLUP", "BayesA", "BayesB", "BayesC", "BL", "BRR", "RKHS")
#   
#   # Error reporting
#   if(is.null(g.in)) { stop("The genotype matrix was not specified.") }
#   if(is.null(y.in)) { stop("The phenotype data was not specified.") }
#   if(nrow(g.in) != nrow(y.in)) { stop("The number of lines in the genotype matrix and in the phenotype data do not match.") }
#   if(!all(models %in% all.models)) {
#     stop("'", paste(models[which(!models %in% all.models)], "' is/are not among the acceptabled models. Please check spelling.", sep = "")) }
#   
#   # Create a matrix of the combinations of trait and model
#   instructions <- expand.grid(models, colnames(y.in))
#   colnames(instructions) <- c("Model", "Trait")
#   # Add other columns, but only include NA as values
#   instructions$r.mean <- NA
#   instructions$r.sd <- NA
#   
#   # Create a function that detects the model and trait and implement cross-validation
#   implement.xval <- function(g.in = g.in, y.in = y.in, trait = trait, model = model, reps = reps, frac.train = frac.train) {
#     
#     # Estimate the number of lines to sample for the training set
#     n.lines <- round(frac.train*nrow(g.in))
#     
#     # Create functions for models
#     cv.rrblup <- function(y.train = y.train, g.train = g.train, model = model) {
#       solve.out <- mixed.solve(y = y.train,  Z = g.train, SE = F, return.Hinv = F)
#       return(solve.out$u)
#     }
#     cv.blr <- function(y.train = y.train, g.train = g.train, model = model) {
#       bayes.out <- BGLR(y = y.train, ETA = list(list(X = g.train, model = model)), verbose = FALSE, nIter = CV.nIter, burnIn = CV.burnIn, saveAt = paste(model, "_", sep = ""))
#       return(bayes.out$ETA[[1]]$b)
#     }
#     
#     # Create a list of the functions so they can be easily referenced
#     cv.func.list <- list(rrBLUP = cv.rrblup, BLR = cv.blr)
#     
#     # Determine the model class (rrBLUP or BLR) for proceeding
#     if (model == "rrBLUP") { model.class <- "rrBLUP" 
#     } else { model.class <- "BLR" }
#     
#     # Notify which trait is being cross-validated
#     writeLines(paste("\nNow conducting fractional cross validation on ", trait, " using ", model, ".", sep = ""))
#     # Start progress bar
#     pb <- txtProgressBar(min = 0, max = reps, style = 3)
#     
#     # Correlation vector
#     corr.vect <- numeric()
# 
#     # For each rep, conduct cross validation
#     for (r in 1:reps) {
#       # Split up the data
#       train <- row.names(g.in)[sample(1:nrow(g.in), n.lines)] # Names of training lines
#       pred <- setdiff(row.names(g.in), train)
#       g.train <- g.in[train,] # Set training genos
#       g.pred <- g.in[pred,] # Set prediction genos
#       y.train <- as.matrix(y.in[train,trait])
#       y.pred <- as.matrix(y.in[pred,trait])
#       
#       mkr.effects <- cv.func.list[model.class][[1]](y.train, g.train, model)
#       
#       # Calculate GEBVs and correlations
#       GEBV <- g.pred %*% mkr.effects
#       corr.vect[r] <- cor(GEBV, y.pred, use = "complete.obs") # correlate GEBVs to actual phenos
#       
#       # Update progress bar
#       setTxtProgressBar(pb = pb, value = r)
#     }
#     
#     # Return the correlation mean and standard deviation
#     results <- cbind(mean(corr.vect), sd(corr.vect))
#     colnames(results) <- c("r.mean", "r.sd")
#     return(results)
#   }
#   
#   # Determine whether to run in parallel or not
#   if(is.null(n.cores)) {
#     
#     # Iterate through the "instructions" data.frame using a for loop
#     for (w in 1:nrow(instructions)) {
#       instructions[w,c(3,4)] <- implement.xval(g.in, y.in, trait = instructions$Trait[w], model = instructions$Model[w], reps, frac.train)
#     
#     }} else {
#       library(parallel) # Load the parallel library
#       # Check if the user is running Windows
#       if (Sys.info()['sysname'] == "Windows") { stop("You are using a Windows system. It is advised to not parallelize R functions in a Windows system.")
#       } else {
#         # Check if the specified number of cores is greater than the cores available
#         if (n.cores > detectCores()) { stop ("The requested number of cores exceeds the number available. Use the 'detectCores' function to determine the number of available cores.")
#         } else {
#           
#           # Break up the instructions data.frame into 'n.cores' equal parts
#           instructions.IT <- split(x = 1:nrow(instructions), f = factor(cut(x = 1:nrow(instructions), breaks = n.cores, labels = FALSE)))
#           
#           # Run cross validation in parallel
#           xval.pieces <- mclapply(X = instructions.IT, FUN = function(x) {
#             instructions.run <- instructions[x,]
#             for (w in 1:nrow(instructions.run)) {
#               instructions.run[w,c(3,4)] <- implement.xval(g.in, y.in, trait = instructions.run$Trait[w], model = instructions.run$Model[w], reps, frac.train)
#             }
#             return(instructions.run)
#           }, mc.cores = n.cores)
#           
#           # Combine the dataframes from parallelization
#           xval.results <- do.call("rbind", xval.pieces)
#         }}}
#   
#   # Compile data into a list
#   xval.results <- list(results = xval.results, models = models, reps = reps, frac.train = frac.train, CV.nIter = CV.nIter, CV.burnIn = CV.burnIn, n.cores = n.cores)
#   return(xval.results)
#   
# } # Close the function


# Define the function for fractional cross-validation (i.e. 60:40, 50:50, etc)
frac.xval <- function(g.in = NULL, # A n x m genotype matrix (coded as -1, 0, 1), with the row.names of g.in equal to the line names
                      y.in = NULL, # A n x t matrix of phentypes, with the row.names of y.in equal to the line names
                      frac.train = 0.70, # The fraction of the lines to use as the training set
                      reps = 250, # The number of iterations of fractional cross-validation
                      model = "rrBLUP"
) {

  # Error reporting
  if(is.null(g.in)) { stop("The genotype matrix was not specified.") }
  if(is.null(y.in)) { stop("The phenotype data was not specified.") }
  if(nrow(g.in) != nrow(y.in)) { stop("The number of lines in the genotype matrix and in the phenotype data do not match.") }
  
  # Create functions for models
  cv.rrblup <- function(y.train = y.train, g.train = g.train) {
    solve.out <- mixed.solve(y = y.train,  Z = g.train, SE = F, return.Hinv = F)
    return(solve.out$u)
  }
  
  # Create a function that detects the model and trait and implement cross-validation
  implement.xval <- function(g.in, y.in, trait, reps, frac.train, model) {
    
    # Estimate the number of lines to sample for the training set
    n.lines <- round(frac.train*nrow(g.in))

    # Notify which trait is being cross-validated
    writeLines(paste("\nNow conducting fractional cross validation on ", trait, " using ", model, ".", sep = ""))
    
    # For each rep, conduct cross validation and output a vector of correlations
    acc.out <- sapply(X = 1:reps, FUN = function(x) {
    
      # Split up the data
      train <- row.names(g.in)[sample(1:nrow(g.in), n.lines)] # Names of training lines
      pred <- setdiff(row.names(g.in), train)
      g.train <- g.in[train,] # Set training genos
      g.pred <- g.in[pred,] # Set prediction genos
      y.train <- as.matrix(y.in[train,trait])
      y.pred <- as.matrix(y.in[pred,trait])
      
      # Calculate marker effects
      u.hat <- cv.rrblup(y.train = y.train, g.train = g.train)
      
      # Calculate GEBVs and correlations
      GEBV <- g.pred %*% u.hat
      acc <- cor(GEBV, y.pred, use = "complete.obs") # correlate GEBVs to actual phenos
      
      return(acc)
    })
    
    # Return the mean and sd of the correlations
    r.mean <- mean(acc.out)
    r.sd <- sd(acc.out)
    return(cbind(trait, model, r.mean, r.sd))
    
  } # Close the function
  
  # Create the output matrix
  results.out <- as.data.frame(matrix(nrow = ncol(y.in), ncol = 4))
  colnames(results.out) <- c("trait", "model", "r.mean", "r.sd")
  
  # Iterate through each trait and perform cross-validation
  for (t in 1:length(colnames(y.in))) {
    trait <- colnames(y.in)[t]
    results.out[t,] <- implement.xval(g.in = g.in, y.in = y.in, trait = trait, reps = reps, frac.train = frac.train, model = model)
  }
  
  # Return a list of results
  cross.val.results <- list(xval.result = results.out, frac.train = frac.train, reps = reps, model = model)
  return(cross.val.results)

} # Close the function




# Define the function for k-fold cross validation (i.e. 10-fold)
k.xval <- function(g.in = NULL, # A n x m genotype matrix (coded as -1, 0, 1), with the row.names of g.in equal to the line names
                   y.in = NULL, # A n x t matrix of phentypes, with the row.names of y.in equal to the line names
                   k.fold = 10, # The number of "folds" to perfom cross-validation with
                   reps = 25, # The number of iterations to divide the lines into k "folds"
                   model = "rrBLUP"
                   ) {
  
  # Error reporting
  if(is.null(g.in)) { stop("The genotype matrix was not specified.") }
  if(is.null(y.in)) { stop("The phenotype data was not specified.") }
  if(nrow(g.in) != nrow(y.in)) { stop("The number of lines in the genotype matrix and in the phenotype data do not match.") }
  
  # Create a function for using rrBLUP
  cv.rrblup <- function(y.train = y.train, g.train = g.train) {
    solve.out <- mixed.solve(y = y.train,  Z = g.train, SE = F, return.Hinv = F)
    return(solve.out$u)
  }
  
  # Create a function that detects the model and trait and implement cross-validation
  implement.xval <- function(g.in, y.in, trait, reps, k.fold, model) {
    
    n.lines <- nrow(g.in)
    
    # Notify which trait is being cross-validated
    writeLines(paste("\nNow conducting fractional cross validation on ", trait, " using ", model, ".", sep = ""))
    
    # sapply function of all reps of k.fold cv for a trait
    acc.out <- sapply(X = 1:reps, FUN = function(x) {
    
      # Design a list of folds
      k.it <- split(x = sample(x = 1:n.lines, size = n.lines), f = factor(cut(1:n.lines, breaks = k.fold)))
      # Apply the function of the list of folds
      k.fold.out <- sapply(X = k.it, FUN = function(x) {
        pred <- row.names(g.in)[x] # Names of prediction set
        train <- setdiff(row.names(g.in), pred) # Names of training set
        g.train <- g.in[train,] # Set training genos
        g.pred <- g.in[pred,] # Set prediction genos
        y.train <- as.vector(y.in[train,trait])
        
        # Marker effects
        u.hat <- cv.rrblup(y.train = y.train, g.train = g.train)
        
        # GEBVs and correlation
        GEBV <- g.pred %*% u.hat
        return(GEBV) # correlate GEBVs to actual phenos
      })
      
      # Create a vector of GEBVs
      GEBV <- do.call("rbind", k.fold.out)
      # Sort it by lines
      GEBV <- GEBV[order(row.names(GEBV))]
      
      # Calculate correlation
      acc <- cor(GEBV, y.in[,trait], use = "complete.obs")
      return(acc)
    })
    
    # Return the mean and sd of the correlations
    r.mean <- mean(acc.out)
    r.sd <- sd(acc.out)
    return(cbind(trait, model, r.mean, r.sd))
    
  } # Close the function
    
  # Create the output matrix
  results.out <- as.data.frame(matrix(nrow = ncol(y.in), ncol = 4))
  colnames(results.out) <- c("trait", "model", "r.mean", "r.sd")
  
  # Iterate through each trait and perform cross-validation
  for (t in 1:length(colnames(y.in))) {
    trait <- colnames(y.in)[t]
    results.out[t,] <- implement.xval(g.in = g.in, y.in = y.in, trait = trait, reps = reps, k.fold = k.fold, model = "rrBLUP")
  }
  
  # Return a list of results
  cross.val.results <- list(xval.result = results.out, folds = k.fold, reps = reps, model = model)
  return(cross.val.results)
} # Close the function
