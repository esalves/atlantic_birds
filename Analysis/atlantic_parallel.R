# atlantic_parallel.R
# ---------------------------------------------------------------------------
# Re-fits the primary wing-length model across a sample of phylogenetic trees
# to propagate phylogenetic uncertainty, then pools across trees.
#
# UPGRADES adopted from the lab's Bergmann/migration manuscript
# (Mizuno, Lundgren, ..., Santos & Nakagawa), see CODE_REVIEW.md §A:
#   (1) PHYLOGENY: trees come from McTavish et al. (2025) "A complete and
#       dynamic tree of birds" via the clootl package, replacing the older
#       Jetz et al. (2012) BirdTree/Hackett backbone (read.nexus of
#       Hackett_trees.nex). clootl is versioned, so record the version used.
#   (2) SPECIES/TREE MATCHING: use prepR4pcm (Nakagawa, Ortega, Mizuno, Santos
#       et al. 2026) to (a) retrieve the tree cloud via its clootl backend and
#       (b) reconcile the trait-data species names to the tree tips through a
#       4-stage cascade (exact -> normalised -> synonym -> fuzzy), with a full
#       audit, so species are not silently dropped. See CODE_REVIEW.md §D.
#   (3) PHYLOGENETIC UNCERTAINTY: sample 50 trees from a pool of 100 and pool
#       parameter estimates ACROSS trees with Rubin's rules
#       (Nakagawa & de Villemereuil 2019) rather than naively concatenating
#       posterior draws with combine_models().
#
# HARMONISED RESPONSE: conc.wing.length (coalesced wing length) for all wing
# models (CODE_REVIEW.md §1.4). MODERN brms SYNTAX: data2 + gr() (§2.1).
#
# Run on a machine with R (>= 4.4), a working Stan toolchain, and enough cores.
# Install: pak::pak("itchyshin/prepR4pcm"); install.packages(c("clootl","phytools")).
# Saves brm0_multiphylo.rda + the Rubin-pooled summary.
# ---------------------------------------------------------------------------

library(brms)
library(ape)
library(MCMCglmm)
library(prepR4pcm)     # species/tree reconciliation + tree retrieval (clootl backend)
library(phytools)      # force.ultrametric()
library(future.apply)
library(posterior)
library(dplyr)

set.seed(20240101)

# --- robust path resolution -------------------------------------------------
# Locate repo files regardless of getwd(): check the working dir, an Analysis/
# subdir, and parent directories (and their Analysis/ subdir) walking upward.
.find_file <- function(fname, subdir = "Analysis") {
  dirs <- getwd(); d <- getwd()
  for (i in 1:8) { d <- dirname(d); dirs <- c(dirs, d) }
  cand <- unique(c(file.path(dirs, fname), file.path(dirs, subdir, fname)))
  hit  <- cand[file.exists(cand)]
  if (!length(hit)) stop("Could not locate '", fname,
      "'. Run from inside the atlantic_birds repo, or setwd() to the Analysis/ folder.")
  normalizePath(hit[1], winslash = "/")
}
ANALYSIS_DIR <- dirname(.find_file("passer90.rda"))  # = .../atlantic_birds/Analysis
apath <- function(...) file.path(ANALYSIS_DIR, ...)  # build all paths from here

load(apath("passer90.rda"))                # cleaned data from atlantic_birds_ms.Rmd
# clootl/eBird uses scientific names WITHOUT underscores; provide a clean column.
passer90$species_name <- gsub("_", " ", passer90$Binomial)

# --- ABT (2018) -> current eBird/Clements names ---------------------------
# clootl errors during retrieval on names absent from the eBird taxonomy, so we
# update genus/epithet changes BEFORE querying. These seven are the genera that
# changed since ABT was compiled (verified against 2024-2025 eBird/Clements +
# SACC). Add to this map if pr_get_tree(..., n_tree = 1)$unmatched flags more.
ebird_synonyms <- c(
  "Antilophia galeata"       = "Chiroxiphia galeata",   # Antilophia merged into Chiroxiphia (2024)
  "Tachyphonus cristatus"    = "Loriotus cristatus",    # moved to Loriotus
  "Pyrrhocoma ruficeps"      = "Thlypopsis pyrrhocoma", # genus AND epithet changed
  "Pyriglena pernambucensis" = "Pyriglena leuconota",   # now a subspecies of leuconota
  "Tangara sayaca"           = "Thraupis sayaca",
  "Tangara cayana"           = "Stilpnia cayana",        # moved to Stilpnia (verify)
  "Tiaris fuliginosus"       = "Asemospiza fuliginosa"   # Tiaris split -> Asemospiza (verify)
)
hit <- passer90$species_name %in% names(ebird_synonyms)
passer90$species_name[hit] <- ebird_synonyms[passer90$species_name[hit]]
spp_data <- unique(passer90$species_name)

# Sanity check before the expensive 100-tree pull (n_tree = 1 reports unmatched
# names instead of erroring); investigate anything it lists.
stopifnot(length(pr_get_tree(spp_data, source = "clootl", n_tree = 1)$unmatched) == 0)

# --- (1) Retrieve the tree cloud with prepR4pcm (clootl backend) ----------
# pr_get_tree() wraps clootl and returns a posterior of trees when n_tree > 1.
# Sampling across trees needs the AvesData repo (the 100-tree dated sample sets):
# download it once (creates <dir>/AvesDataLite-main) and point clootl at it.
# Ensure the AvesData repo is available. clootl records its location in the
# AVESDATA_PATH environment variable (written to .Renviron by get_avesdata_repo),
# so only download when that is unset/missing — this avoids re-downloading when
# the repo already lives elsewhere (e.g. you ran get_avesdata_repo(path=".") from
# the repo root). The download is the dominant cost, so this guard matters.
if (!nzchar(Sys.getenv("AVESDATA_PATH")) || !dir.exists(Sys.getenv("AVESDATA_PATH"))) {
  clootl::get_avesdata_repo(path = ANALYSIS_DIR)
}
# cache = TRUE stores the result on disk and reuses it on identical calls, so
# re-sourcing this script does not re-query clootl every time.
got    <- pr_get_tree(spp_data, source = "clootl", n_tree = 100, cache = TRUE)
trees  <- got$tree                              # multiPhylo (100 dated trees)
clootl_version <- pr_cite_tree(got, format = "text")  # provenance for the methods

# --- (2) Reconcile data species names to the tree tips --------------------
# Four-stage cascade (exact -> normalised -> synonym -> fuzzy), fully audited,
# so no species is silently dropped. Review the result before applying.
rec <- reconcile_tree(
  x         = passer90,
  tree      = trees[[1]],                       # taxa are identical across the cloud
  x_species = "species_name",
  fuzzy     = TRUE,
  resolve   = "flag"
)
print(reconcile_summary(rec))                   # inspect coverage + each match stage
# reconcile_report(rec)                         # full provenance, if you want it
# Genus/epithet changes were already fixed upfront via `ebird_synonyms`, so this
# pass should resolve all species. For any residual mismatch (typo/synonym not in
# that map), fix here, e.g.:
#   rec <- reconcile_override(rec, name_x = "Xxx yyy", name_y = "<eBird tip>")
#   # or use the bundled crosswalk: reconcile_crosswalk(crosswalk_birdlife_birdtree)
# For a species genuinely absent from the tree, reconcile_augment() can graft it
# as sister to a congener — then fit models WITH and WITHOUT the graft and report
# whether conclusions change (see CODE_REVIEW.md §D).

# Aligned data + pruned tree, with the data's species column harmonised to the
# tip labels (identical species sets — the precondition for any PCM).
aligned <- reconcile_apply(rec, data = passer90, tree = trees[[1]],
                           species_col = "species_name", drop_unresolved = TRUE)
dat <- aligned$data
dat$species_name <- gsub(" ", "_", dat$species_name)  # underscore form for brms grouping
dat$spp <- dat$species_name                            # both grouping terms use matched names

# Prune ALL 100 trees to the resolved species set, then subsample 50.
# reconcile_apply() relabels the aligned tree's tips with spaces, but the clootl
# cloud uses underscores; normalise both to underscores and keep the
# INTERSECTION so keep.tip() never receives an unmatched label.
norm_us   <- function(x) gsub(" ", "_", x)
keep_tips <- norm_us(aligned$tree$tip.label)
trees_pruned <- lapply(trees, function(t) {
  t$tip.label <- norm_us(t$tip.label)
  ape::keep.tip(t, intersect(keep_tips, t$tip.label))
})
class(trees_pruned) <- "multiPhylo"
tree_samp <- sample(trees_pruned, 50)           # seed set above

# --- helper: per-tree phylogenetic covariance matrix ----------------------
make_A <- function(tree) {
  tree$tip.label <- gsub(" ", "_", tree$tip.label)  # match dat$species_name format
  # clootl trees are dated but pruning + floating-point leaves them very slightly
  # non-ultrametric, which inverseA(scale = TRUE) refuses. Coerce to ultrametric
  # (minimal NNLS adjustment) first.
  if (!ape::is.ultrametric(tree)) tree <- phytools::force.ultrametric(tree, method = "nnls")
  inv <- inverseA(tree, nodes = "TIPS", scale = TRUE)
  A   <- solve(inv$Ainv); rownames(A) <- rownames(inv$Ainv); A
}

priors <- c(prior(normal(71, 15), class = Intercept),
            prior(normal(0, 10),  class = b),
            prior(cauchy(0, 1),   class = sd),
            prior(cauchy(0, 1),   class = sigma))

fit_one <- function(tree) {
  A <- make_A(tree)
  brm(conc.wing.length ~ 1 + Sex + scaled_yr + scaled_lat +
        (1 + scaled_yr || spp) + (1 | gr(species_name, cov = A)),
      data = dat, data2 = list(A = A),
      family = gaussian(), prior = priors,
      iter = 1500, warmup = 1000, chains = 1, cores = 1,
      control = list(max_treedepth = 12, adapt_delta = 0.95),
      seed = 20240101)
}

plan(multisession)
fits <- future_lapply(tree_samp, fit_one, future.seed = TRUE)

# --- (2) Rubin's rules pooling across trees -------------------------------
# For each fixed effect: pooled estimate = mean of per-tree posterior means;
# total variance = within-tree var (mean of per-tree posterior variances)
# + between-tree var (variance of per-tree means) inflated by (1 + 1/m).
# (Nakagawa & de Villemereuil 2019, Syst. Biol. 68:632-641.)
pool_rubin <- function(fits, pars = c("b_Intercept","b_SexMale","b_scaled_yr","b_scaled_lat")) {
  m <- length(fits)
  draws <- lapply(fits, function(f) as_draws_df(f)[, pars, drop = FALSE])
  means <- t(sapply(draws, function(d) sapply(d, mean)))
  vars  <- t(sapply(draws, function(d) sapply(d, var)))
  qbar  <- colMeans(means)                 # pooled point estimate
  ubar  <- colMeans(vars)                  # within-imputation variance
  b     <- apply(means, 2, var)            # between-tree variance
  tot   <- ubar + (1 + 1/m) * b            # total variance
  se    <- sqrt(tot)
  data.frame(par = pars, estimate = qbar, se = se,
             lower = qbar - 1.96*se, upper = qbar + 1.96*se)
}

rubin_summary <- pool_rubin(fits)
print(rubin_summary)

save(fits, rubin_summary, clootl_version, file = apath("brm0_multiphylo.rda"))

# NOTE: Pagel's-lambda-style variance partitioning and parametric-bootstrap
# CIs for the variance components (as in the Mizuno et al. paper) can be added
# with a phylolm fit on a consensus tree rescaled to unit tip height; see
# CODE_REVIEW.md §A for the rationale and the brms vs phylolm trade-off.
