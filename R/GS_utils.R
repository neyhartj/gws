## GS_utils.R
# An R script containing miscellaneous functions for genomic prediction
# Author: Jeff Neyhart
# Date: January 31, 2016

## vp.predict
# A function to use phenotypic and genotypic data from a training population, along with
# genotypic data of a validation population, to predict the phenotypes of the validation
# population using RR-BLUP. Also uses known phenotype data of a validation population to 
# estimate the prediction accuracy. 
# This function is not error-optimized, so use at own risk
## Variables to note
# nt = number of training genotypes
# nv = number of validation genotypes
# m = number of markers
# p = number of traits

# Load libraries
library(boot)
library(EMMREML)

vp.predict <- function(y.train, # nt x 1 vector of training phenotypes
                       g.train, # nt x m matrix of training marker information
                       y.pred, # nv x 1 vector of validation phenotypes
                       g.pred, # nv x m matrix of validation marker information
                       bootreps = 1000, # Number of replicates of bootstrapping for SE estimation
                       PEV = TRUE # Logical scalar to return the prediction error variance
                       ) {

  # Define mixed model variables
  y <- y.train
  Z <- g.train
  
  # Solve the mixed model
  solve.out <- mixed.solve(y = y, Z = Z, method = "REML", return.Hinv = FALSE)
  # Extract marker effects
  marker.effects <- solve.out$u
  # Calculate GEBVs
  GEBVs <- g.pred %*% marker.effects
  
  # Calculate a bootstrapped correlation
  bootdata <- cbind(GEBVs, y.pred)
  bootcor <- function(d, i)  {
    d2 <- d[i,]
    return(cor(d2[,1], d2[,2]))
  }
  boot.results <- boot(bootdata, bootcor, R = bootreps)
  
  # Extract data
  r.pred <- boot.results$t0
  r.pred.sd <- sd(boot.results$t)
  
  return(cbind(r.pred, r.pred.sd))
} # Close the function