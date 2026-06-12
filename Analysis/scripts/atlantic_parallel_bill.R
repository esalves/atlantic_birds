# atlantic_parallel_bill.R
# ---------------------------------------------------------------------------
# Bill-width counterpart of atlantic_parallel.R. Fits the bill-width model
# across the SAME pipeline as the wing model — McTavish/clootl tree cloud,
# prepR4pcm name reconciliation, 50 trees, Rubin's-rules pooling — so the bill
# result is reported on the same updated phylogeny as wing length
# (CODE_REVIEW.md §A/§D; resolves the bill-width [AUTHOR ACTION] in index.qmd).
#
# Differences from the wing model (matching the manuscript's design):
#   * response  = Bill_width.mm.   (a single measured column; NA rows dropped)
#   * NO Sex term
#   * intercept prior Normal(9, 2.5)  (bill width is ~9 mm, not ~71 mm)
#
# Retrieval is cached (cache = TRUE), so this reuses the cloud already pulled by
# atlantic_parallel.R; only the 50 bill fits cost time (fewer obs + one fewer
# parameter than wing, so faster). Saves brm_bill_multiphylo.rda.
#
# Install: pak::pak("itchyshin/prepR4pcm"); install.packages(c("clootl","phytools")).
# ---------------------------------------------------------------------------

library(brms)
library(ape)
library(MCMCglmm)
library(prepR4pcm)
library(phytools)
library(future.apply)
library(posterior)
library(dplyr)

set.seed(20240101)

# --- robust path resolution (repo layout: Analysis/{scripts,data,output,figures}) ---
.find_analysis_dir <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:10) {
    if (dir.exists(file.path(d, "data", "derived")) && dir.exists(file.path(d, "scripts")))
      return(d)
    if (dir.exists(file.path(d, "Analysis", "data", "derived")))
      return(normalizePath(file.path(d, "Analysis"), winslash = "/"))
    parent <- dirname(d); if (identical(parent, d)) break; d <- parent
  }
  stop("Could not locate Analysis/ (need Analysis/data/derived and Analysis/scripts). ",
       "Run from inside the atlantic_birds repo.")
}
ANALYSIS_DIR <- .find_analysis_dir()
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
fig_path     <- function(...) file.path(ANALYSIS_DIR, "figures", ...)
script_path  <- function(...) file.path(ANALYSIS_DIR, "scripts", ...)

load(derived_path("passer90.rda"))
passer90$species_name <- gsub("_", " ", passer90$Binomial)

# ABT (2018) -> current eBird/Clements names (see atlantic_parallel.R).
ebird_synonyms <- c(
  "Antilophia galeata"       = "Chiroxiphia galeata",
  "Tachyphonus cristatus"    = "Loriotus cristatus",
  "Pyrrhocoma ruficeps"      = "Thlypopsis pyrrhocoma",
  "Pyriglena pernambucensis" = "Pyriglena leuconota",
  "Tangara sayaca"           = "Thraupis sayaca",
  "Tangara cayana"           = "Stilpnia cayana",
  "Tangara palmarum"         = "Thraupis palmarum",
  "Tangara peruviana"        = "Stilpnia peruviana",
  "Dixiphia pipra"           = "Pseudopipra pipra",
  "Tiaris fuliginosus"       = "Asemospiza fuliginosa"
)
hit <- passer90$species_name %in% names(ebird_synonyms)
passer90$species_name[hit] <- ebird_synonyms[passer90$species_name[hit]]
spp_data <- unique(passer90$species_name)

# --- tree retrieval + reconciliation (identical to the wing pipeline) ------
# Only download the AvesData repo if clootl has no recorded location for it
# (AVESDATA_PATH); avoids a redundant re-download when it already exists elsewhere.
if (!nzchar(Sys.getenv("AVESDATA_PATH")) || !dir.exists(Sys.getenv("AVESDATA_PATH"))) {
  clootl::get_avesdata_repo(path = raw_path())
}
.check <- pr_get_tree(spp_data, source = "clootl", n_tree = 1)
if (length(.check$unmatched) > 0) {
  message("Removing ", length(.check$unmatched), " species absent from eBird taxonomy: ",
          paste(.check$unmatched, collapse = ", "))
  spp_data <- setdiff(spp_data, .check$unmatched)
}
got   <- pr_get_tree(spp_data, source = "clootl", n_tree = 100, cache = TRUE)
trees <- got$tree

rec     <- reconcile_tree(passer90, trees[[1]], x_species = "species_name",
                          fuzzy = TRUE, resolve = "flag")
print(reconcile_summary(rec))
aligned <- reconcile_apply(rec, data = passer90, tree = trees[[1]],
                           species_col = "species_name", drop_unresolved = TRUE)
dat <- aligned$data
dat$species_name <- gsub(" ", "_", dat$species_name)
dat$spp <- dat$species_name

norm_us   <- function(x) gsub(" ", "_", x)
keep_tips <- norm_us(aligned$tree$tip.label)
trees_pruned <- lapply(trees, function(t) {
  t$tip.label <- norm_us(t$tip.label)
  ape::keep.tip(t, intersect(keep_tips, t$tip.label))
})
class(trees_pruned) <- "multiPhylo"
set.seed(20240101)
tree_samp <- sample(trees_pruned, 50)

make_A <- function(tree) {
  tree$tip.label <- gsub(" ", "_", tree$tip.label)
  if (!ape::is.ultrametric(tree)) tree <- phytools::force.ultrametric(tree, method = "nnls")
  inv <- inverseA(tree, nodes = "TIPS", scale = TRUE)
  A   <- solve(inv$Ainv); rownames(A) <- rownames(inv$Ainv); A
}

# --- bill-width model ------------------------------------------------------
priors_bill <- c(prior(normal(9, 2.5), class = Intercept),  # bill width ~ 9 mm
                 prior(normal(0, 10),  class = b),
                 prior(cauchy(0, 1),   class = sd),
                 prior(cauchy(0, 1),   class = sigma))

fit_bill <- function(tree) {
  A <- make_A(tree)
  brm(Bill_width.mm. ~ 1 + scaled_yr + scaled_lat +
        (1 + scaled_yr || spp) + (1 | gr(species_name, cov = A)),
      data = dat, data2 = list(A = A),
      family = gaussian(), prior = priors_bill,
      iter = 1500, warmup = 1000, chains = 1, cores = 1,
      control = list(max_treedepth = 12, adapt_delta = 0.95),
      seed = 20240101)
}

plan(multisession)
bill_fits <- future_lapply(tree_samp, fit_bill, future.seed = TRUE)

# --- Rubin's rules pooling across trees (no SexMale term here) --------------
pool_rubin <- function(fits, pars = c("b_Intercept", "b_scaled_yr", "b_scaled_lat")) {
  m <- length(fits)
  draws <- lapply(fits, function(f) as_draws_df(f)[, pars, drop = FALSE])
  means <- t(sapply(draws, function(d) sapply(d, mean)))
  vars  <- t(sapply(draws, function(d) sapply(d, var)))
  qbar  <- colMeans(means)
  ubar  <- colMeans(vars)
  b     <- apply(means, 2, var)
  tot   <- ubar + (1 + 1/m) * b
  se    <- sqrt(tot)
  data.frame(par = pars, estimate = qbar, se = se,
             lower = qbar - 1.96*se, upper = qbar + 1.96*se)
}
bill_rubin_summary <- pool_rubin(bill_fits)
print(bill_rubin_summary)

# phylogenetic signal (proportion of variance), pooled across the 50 fits
bill_sig <- sapply(bill_fits, function(f) {
  hypothesis(f,
    "sd_species_name__Intercept^2 / (sd_species_name__Intercept^2 + sd_spp__Intercept^2 + sd_spp__scaled_yr^2 + sigma^2) = 0",
    class = NULL)$hypothesis$Estimate
})
print(round(c(mean = mean(bill_sig), lwr = quantile(bill_sig, .025), upr = quantile(bill_sig, .975)), 3))

save(bill_fits, bill_rubin_summary, bill_sig, file = out_path("models", "brm_bill_multiphylo.rda"))
