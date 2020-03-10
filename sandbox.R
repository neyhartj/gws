# # Sandbox for gws
#
#
# # Compare pop_predict_quick with pop.predict
#
# library(tidyverse)
# library(gws)
# library(PopVar)
#
# # Load data
# data("genos")
# data("phenos")
# data("map")
#
# # Create a crossing block
# crossing.table <- combn(x = row.names(genos), m = 2) %>%
#   t() %>%
#   as.data.frame() %>%
#   structure(names = c("parent1", "parent2")) %>%
#   sample_n(100)
#
# # All traits
# y_use <- phenos
#
# traits <- colnames(y_use)[-1]
#
# ## Round the genotypes
# genos1 <- ifelse(genos < 0, -1, 1)
# # Convert for use in PopVar
# genos_popvar <- as.data.frame(cbind( c("", row.names(genos1)), rbind(colnames(genos1), genos1)) )
#
#
# ## rename for function inputs
# G.in <- genos1
# y.in <- y_use
# map.in <- map
#
# #
#
# pp_quick_out <- pop_predict_quick(G.in = G.in, y.in = y.in, map.in = map.in, crossing.table = crossing.table)
#
#
#
# ## Alternate for PopVar
# G.in <- genos_popvar
#
# pp_out <- pop.predict(G.in = G.in, y.in = y.in, map.in = map.in, crossing.table = crossing.table,
#                       min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, impute = "pass", nInd = 1000, nSim = 1,
#                       nCV.iter = 1, models = "rrBLUP", remove.dups = FALSE)
#
#
# ## Predictions using pop.predict2
# pp2_out <- pop.predict2(G.in = G.in, y.in = y.in, map.in = map.in, crossing.table = crossing.table, self.gen = Inf)
#
#
#
# # Tidy
# if (length(traits) > 1) {
#
#   pp_out_tidy <- pp_out$predictions %>%
#     do.call("rbind", .) %>%
#     as.data.frame() %>%
#     mutate_all(unlist) %>%
#     mutate(trait = rep(traits, each = nrow(.) / length(traits)))
#
# } else {
#
#   pp_out_tidy <- pp_out$predictions %>%
#     as.data.frame() %>%
#     mutate_all(unlist) %>%
#     mutate(trait = rep(traits, each = nrow(.) / length(traits)))
#
# }
#
#
# # Test correlations
# pp_quick_out_FHB <- pp_quick_out %>%
#   filter(trait == "FHB")
#
# pp_out_tidy_FHB <- pp_out_tidy %>%
#   filter(trait == "FHB")
#
# pp2_out_tidy_FHB <- pp2_out %>%
#   filter(trait == "FHB")
#
#
# # Comparison of predicted means
# c(acc_pred_mu = cor(pp_quick_out_FHB$pred_mu, pp_out_tidy_FHB$pred.mu),
#   acc_pred_varG = cor(pp_quick_out_FHB$pred_varG, pp_out_tidy_FHB$pred.varG),
#   acc_pred_mu_sp = cor(pp_quick_out_FHB$pred_mu_sp_high, pp_out_tidy_FHB$mu.sp_high),
#   acc_pred_corG = cor(pp_quick_out_FHB$DON, pp_out_tidy_FHB$`cor_w/_DON`))
#
# c(acc_pred_mu = cor(pp2_out_tidy_FHB$pred_mu, pp_out_tidy_FHB$pred.mu),
#   acc_pred_varG = cor(pp2_out_tidy_FHB$pred_varG, pp_out_tidy_FHB$pred.varG),
#   acc_pred_mu_sp = cor(pp2_out_tidy_FHB$pred_musp_high, pp_out_tidy_FHB$mu.sp_high),
#   acc_pred_corG = cor(pp2_out_tidy_FHB$cor_w_DON, pp_out_tidy_FHB$`cor_w/_DON`))
#
#
#
#
#
#
# ## Plot
# par(mfrow = c(2,2))
# plot(pp_quick_out_FHB$pred_mu, pp_out_tidy_FHB$pred.mu, main = "Comparison of predicted family means",
#      xlab = "Deterministic", ylab = "PopVar")
# plot(pp_quick_out_FHB$pred_varG, pp_out_tidy_FHB$pred.varG, main = "Comparison of predicted genetic variance",
#      xlab = "Deterministic", ylab = "PopVar")
# plot(pp_quick_out_FHB$pred_mu_sp_low, pp_out_tidy_FHB$mu.sp_low, main = "Comparison of predicted superior progeny means",
#      xlab = "Deterministic", ylab = "PopVar")
# plot(pp_quick_out_FHB$Yield, pp_out_tidy_FHB$`cor_w/_Yield`,
#      main = "Comparison of predicted genetic correlations between\nFHB Severity and Yield",
#      xlab = "Deterministic", ylab = "PopVar")
# plot(pp2_out_tidy_FHB$cor_w_Yield, pp_out_tidy_FHB$`cor_w/_Yield`,
#      main = "Comparison of predicted genetic correlations between\nFHB Severity and Yield",
#      xlab = "Deterministic", ylab = "PopVar")
#
#
#
#
#
#
#
# ### Debug pop.predict2 from Aaron
# library(gws)
# library(tidyverse)
#
#
# ## Read in the files
# blues5 <- read.csv(file = "C:/Users/jln54/Downloads/blues5.csv")
# crossTable <- read.csv(file = "C:/Users/jln54/Downloads/crossTable.csv")
# genoAll3 <- read.csv(file = "C:/Users/jln54/Downloads/genoAll3_1.csv", header = FALSE)
# genoForMap2 <- read.csv(file = "C:/Users/jln54/Downloads/genoForMap2.csv", na.strings = c("NA", "#N/A"))
#
# pop_predict_out <- pop.predict2(G.in = genoAll3, y.in = blues5, map.in = genoForMap2, crossing.table = crossTable)
#
# pop_predict_out <- pop.predict2(G.in = genoAll3[,-1], y.in = blues5[,-1], map.in = genoForMap2[,-1], crossing.table = crossTable[,-1])
#
#
#
# ## Adjustments to ensure operations
# map.in <- subset(genoForMap2[,-1], !is.na(pos))
#
# genoAll3_1 <- as.data.frame(genoAll3[,-1])
# G.in <- rbind(colnames(genoAll3_1), genoAll3_1)
# G.in1 <- G.in[,c(1, which(unlist(G.in[1, , drop = TRUE]) %in% map.in$mrk))]
#
#
# y.in <- blues5[,-1]
#
#
#
# crossing.table <- crossTable[,-1]
#
# ## Test pop.predict2
# soy_pop_predict <- pop.predict2(G.in = G.in, y.in = y.in, map.in = map.in, crossing.table = crossing.table)





## Compare deterministic with simulation ##

## MultiParentPopVar
##
## Simulations to compare deterministic formula with stochaistic simulations
##
##

# Libraries
library(tidyverse)
library(pbsim)
library(pbsimData)
library(gws)
library(qtl)
library(mpMap2)

# Project repositories
proj_dir <- "C:/GoogleDrive/Projects/SideProjects/MultiParentPopVar/"
data_dir <- file.path(proj_dir, "Data")


## Simulations parameters
sim_pop_size <- 100 # Population size
nIter <- 25 # simulation reps
nCrosses <- 200 # Number of crosses to simulate
n_mar_per_chrom <- 100 # Markers per chromosome to sample

# Create a map
map_sim <- s2_snp_info %>%
  split(.$chrom) %>%
  map(~setNames(.$cM_pos, .$rs)) %>%
  map(~structure(., class = "A"))  %>%
  structure(., class = "map") %>%
  # Jitter
  qtl::jittermap(.) %>%
  `names<-`(., seq_along(.))

## Sample 10 markers per chromosome
set.seed(1730)

{
  # map_sim_use <- structure(map(map_sim, ~sort(sample(x = ., size = n_mar_per_chrom))), class = "map")
  map_sim_use <- map_sim



  ## Marker genotypes to use as base population
  ## Subset markers sampled in the map
  base_pop_geno <- s2_cap_genos[,unlist(map(map_sim_use, names)),drop = FALSE]

  # Name of these potential parents - sample 10
  parents <- sort(sample(row.names(base_pop_geno), 150))
  parent_geno <- base_pop_geno[parents,,drop = FALSE]

  # Create a crossing block with random 4-way crosses
  # Create the crossing block of parents (in (AxB)x(CxD) order)
  cb_4way <- as.data.frame(t(combn(x = sample(parents, 10), m = 4)), stringsAsFactors = FALSE) %>%
    rename_all(~paste0("parent", 1:4)) %>%
    sample_n(tbl = ., size = nCrosses)

}


## Deterministic predictions ##

# Convert the map for deterministic predictions
map_in_use1 <- map2table(map_sim_use) %>%
  rownames_to_column("marker")

# Prepare genotypes
geno_use <- parent_geno - 1
G.in <- as.data.frame(cbind(row.names(geno_use), geno_use), stringsAsFactors = FALSE)
G.in <- rbind(colnames(G.in), G.in)

# Create random marker effects
marker_u <- matrix(data = rnorm(ncol(geno_use)), ncol = 1, dimnames = list(colnames(geno_use), NULL))

# Genotypic values
g_val <- geno_use %*% marker_u
# Random deviates
p_val <- g_val + rnorm(n = nrow(g_val), mean = 0, sd = sqrt(var(g_val) / 2))

# Create random phenotype df
y.in <- data.frame(line_name = parents, y = p_val, stringsAsFactors = FALSE, row.names = NULL)


# Debugging
map.in = map_in_use1
crossing.table = cb_4way
# crossing.table = cb_4way[,1:2]
self.gen = 10
tail.p = 0.1
DH = FALSE
model = "rrBLUP"



# Feed into mpppop_predict - 4-way crosses
varG_4way_mppop_predict <- mppop_predict(G.in = G.in, y.in = y.in, map.in = map_in_use1,
                                         crossing.table = cb_4way, self.gen = self.gen)



## Predict bi-parental crosses using mppop_predict
cb_2way <- cb_4way[,1:2]
varG_2way_mppop_predict <- mppop_predict(G.in = G.in, y.in = y.in, map.in = map_in_use1,
                                         crossing.table = cb_2way, self.gen = self.gen)

# Using pop.predict
varG_2way_pop_predict <- pop.predict2(G.in = G.in, y.in = y.in, map.in = map_in_use1,
                                      crossing.table = cb_2way, self.gen = self.gen)


## Compare
plot(varG_2way_mppop_predict$pred_mu, varG_2way_pop_predict$pred_mu, main = "pred mu")
plot(varG_2way_mppop_predict$pred_varG, varG_2way_pop_predict$pred_varG, main = "pred varG")



## Simulate bi-parental populations to compare

library(rrBLUP)

# Predict marker effects
mf <- model.frame(y ~ line_name, y.in)
y <- model.response(mf)
Z <- geno_use

fit <- mixed.solve(y = y, Z = Z)

marker_effects_fit <- as.matrix(fit$u)
pgv <-( Z %*% marker_effects_fit) + c(fit$beta)

# Data.frame for comparisons
varG_2way_sim <- bind_cols(
  select(varG_2way_mppop_predict, parent1, parent2 = parent3, pred_varG_mppop = pred_varG),
  select(varG_2way_pop_predict, pred_varG_pop = pred_varG)
) %>% mutate(pred_varG_sim = NA)

varG_4way_sim <- varG_4way_mppop_predict %>%
  select(contains("parent"), pred_varG) %>%
  mutate(pred_varG_sim = NA)


# Pedigree
pedigree_use <- mpMap2::fourParentPedigreeSingleFunnel(initialPopulationSize = sim_pop_size,
                                                       selfingGenerations = self.gen, nSeeds = 1,
                                                       intercrossingGenerations = 0)


# pedigree_use <- mpMap2::rilPedigree(populationSize = sim_pop_size, selfingGenerations = self.gen)



for (j in seq_len(nrow(crossing.table))) {

  # Vector of parent names
  pars <- as.character(crossing.table[j,])

  # Calculate the mean PGV
  pred_mu <- mean(pgv[pars,])

  ## Simulate a population 25 times
  sim_out <- replicate(n = nIter, {

    sim_pop <- mpMap2::simulateMPCross(map = map_sim_use, pedigree = pedigree_use, mapFunction = "haldane")

    # Extract marker data of the progeny
    prog_geno <- sim_pop@geneticData@.Data[[1]]@finals

    # Configure matrix of parent genotypes
    par_geno <- parent_geno[pars,,drop = FALSE] - 1

    ## Overlay parent genotypes on progeny
    prog_geno1 <- apply(X = rbind(par_geno, prog_geno), MARGIN = 2, FUN = function(snp) {
      snp[seq_len(nrow(par_geno))][snp[-seq_len(nrow(par_geno))]] })

    # NAs are hets
    prog_geno1[is.na(prog_geno1)] <- 0
    row.names(prog_geno1) <- row.names(prog_geno)

    # Predict PGVs
    c((prog_geno1 %*% marker_effects_fit) + c(fit$beta))

  })

  ## Summarize mean and variances
  pred_varG_sim <- mean(apply(X = sim_out, MARGIN = 2, FUN = var))

  # Report variance
  # varG_2way_sim$pred_varG_sim[j] <- pred_varG_sim
  varG_4way_sim$pred_varG_sim[j] <- pred_varG_sim

}


## Plot
varG_2way_sim



cor(varG_4way_sim$pred_varG, varG_4way_sim$pred_varG_sim)



































