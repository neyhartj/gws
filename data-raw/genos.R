## code to prepare `genos` dataset goes here

# Get data from PopVar
library(PopVar)

# Recode the marker genotypes
M <- as.matrix(sapply(X = G.in_ex[-1, -1], FUN = function(x) as.numeric(as.character(x))))
dimnames(M) <- list(as.character(G.in_ex[,1][-1]), as.character(unlist(G.in_ex[1,][-1])))

# Set 0 to NA; set to mean
M[M == 0] <- NA
genos <- apply(X = M, MARGIN = 2, FUN = function(snp) {
  snp_mean <- mean(snp, na.rm = TRUE)
  snp[is.na(snp)] <- ifelse(snp_mean > 0, 1, -1)
  snp
})


usethis::use_data(genos, overwrite = TRUE)
