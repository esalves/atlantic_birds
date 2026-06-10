# atlantic_parallel.R
# ---------------------------------------------------------------------------
# Re-fits the primary wing-length model across a sample of phylogenetic trees
# to propagate phylogenetic uncertainty, then combines the fits.
#
# HARMONISED RESPONSE: uses conc.wing.length (coalesced right/left/unspecified
# wing length) for ALL wing models, matching the primary single-tree model
# brm0. (The earlier server run used Wing_length_right.mm. only, which roughly
# halved the sample; see CODE_REVIEW.md item 1.4.)
#
# MODERN brms SYNTAX: phylogenetic covariance is passed via data2 + gr(),
# replacing the removed cov_ranef argument (CODE_REVIEW.md item 2.1).
#
# Run on a machine with R, a current CMDSTAN/rstan toolchain, the BirdTree
# phylogeny (Hackett_trees.nex), and enough cores. Saves brm0_multiphylo.rda.
# ---------------------------------------------------------------------------

library(brms)
library(ape)
library(geiger)
library(MCMCglmm)
library(future.apply)

set.seed(20240101)

load("passer90.rda")        # cleaned data (produced by atlantic_birds_ms.Rmd)
trees <- read.nexus("Hackett_trees.nex")

# species present in the data
spp_data <- unique(passer90$Binomial)

# build a per-tree phylogenetic covariance matrix
make_A <- function(tree) {
  drop <- setdiff(tree$tip.label, spp_data)
  tr   <- drop.tip(tree, drop)
  inv  <- inverseA(tr, nodes = "TIPS", scale = TRUE)
  A    <- solve(inv$Ainv)
  rownames(A) <- rownames(inv$Ainv)
  A
}

n_trees <- 50
tree_samp <- sample(trees, size = n_trees)

# weakly-informative priors (as in the single-tree primary model)
priors <- c(prior(normal(71, 15), class = Intercept),
            prior(normal(0, 10),  class = b),
            prior(cauchy(0, 1),   class = sd),
            prior(cauchy(0, 1),   class = sigma))

fit_one <- function(tree) {
  A <- make_A(tree)
  brm(conc.wing.length ~ 1 + Sex + scaled_yr + scaled_lat +
        (1 + scaled_yr || spp) + (1 | gr(Binomial, cov = A)),
      data   = passer90,
      data2  = list(A = A),
      family = gaussian(),
      prior  = priors,
      iter = 1500, warmup = 1000, chains = 1, cores = 1,
      control = list(max_treedepth = 12, adapt_delta = 0.95),
      seed = 20240101)
}

plan(multisession)   # adjust to your core count
brm0 <- future_lapply(tree_samp, fit_one, future.seed = TRUE)

# combine across trees (each is one chain) for downstream summaries
brm0_multiphylo <- combine_models(mlist = brm0, check_data = FALSE)

save(brm0, brm0_multiphylo, file = "brm0_multiphylo.rda")
print(summary(brm0_multiphylo))
