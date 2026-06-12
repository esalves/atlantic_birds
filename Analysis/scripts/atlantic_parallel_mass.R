# atlantic_parallel_mass.R
# ---------------------------------------------------------------------------
# Body-mass counterpart of atlantic_parallel.R. Asks the same question for
# log(body mass) that the main analysis asks for wing length: does body size
# decline over 1990–2018, on the same updated phylogeny and pooling scheme?
#
# Mirrors atlantic_parallel.R EXACTLY (McTavish/clootl tree cloud, prepR4pcm
# name reconciliation, 50 trees, Rubin's-rules pooling, Sex + scaled_yr +
# scaled_lat fixed effects, uncorrelated species intercept+slope, phylogenetic
# covariance). The ONLY differences:
#   * response   = log(Body_mass.g.)      [wing used conc.wing.length]
#   * intercept prior Normal(3, 1.5)      [log grams ~ 3; wing used Normal(71,15)]
#
# WHY a separate sample build: passer90.rda drops Body_mass.g. in its column
# select, so we re-derive the analytical frame from the raw ABT csv using the
# IDENTICAL filter/scaling order as atlantic_birds_ms.Rmd (Year>=1990, Adult,
# Sex!=Unknown, Passeriformes, AF 20-km buffer; scale year/lat; then species
# n>=30 records and span>=5 yr) and additionally keep body mass. This reproduces
# passer90's 89-species / 15,332-record sample (verified) plus the mass column;
# the model then fits on the non-NA mass records, exactly as the wing model
# fits on non-NA wing records. Saves brm_mass_multiphylo.rda + a results summary.
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
library(readr)

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

# --- build the analytical frame WITH body mass (mirrors the Rmd's filter) ---
# Same operations and order as atlantic_birds_ms.Rmd's `newdataframes` chunk, so
# scaled_yr / scaled_lat and the species set match passer90 exactly; we just keep
# Body_mass.g. (which passer90's select() drops) and add ln_body_mass.
birds <- read_csv(raw_path("ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv"),
                  guess_max = 70000, show_col_types = FALSE)
passer90 <- birds %>%
  filter(Year >= 1990, Age == "Adult", Sex != "Unknown",
         Order == "Passeriformes",
         AtlanticForests_20km_Buffer == "inside the 20 km polygon") %>%
  mutate(
    conc.wing.length = coalesce(Wing_length_right.mm., Wing_length_left.mm., Wing_length.mm.),
    ln_body_mass     = log(Body_mass.g.),
    scaled_lat       = as.numeric(scale(Latitude_decimal_degrees * -1)),
    scaled_yr        = as.numeric(scale(Year)),
    Binomial         = gsub(" ", "_", Binomial),
    spp              = Binomial
  ) %>%
  group_by(Binomial) %>%
  filter(n() >= 30) %>%
  mutate(range = max(Year) - min(Year)) %>%
  filter(range >= 5) %>%
  ungroup() %>%
  as.data.frame()

# Sanity: this must reproduce the main sample (89 species / 15,332 records).
cat(sprintf("Sample: %d records, %d species (expect 15332 / 89).\n",
            nrow(passer90), length(unique(passer90$Binomial))))
cat(sprintf("Body-mass records (non-NA): %d (%.0f%%); log-mass mean %.2f, sd %.2f\n",
            sum(!is.na(passer90$ln_body_mass)),
            100 * mean(!is.na(passer90$ln_body_mass)),
            mean(passer90$ln_body_mass, na.rm = TRUE),
            sd(passer90$ln_body_mass,   na.rm = TRUE)))
save(passer90, file = derived_path("passer90_mass.rda"))

# clootl/eBird uses scientific names WITHOUT underscores; provide a clean column.
passer90$species_name <- gsub("_", " ", passer90$Binomial)

# --- ABT (2018) -> current eBird/Clements names (identical to atlantic_parallel.R) ---
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

# --- tree retrieval + reconciliation (identical pipeline to the wing model) --
if (!nzchar(Sys.getenv("AVESDATA_PATH")) || !dir.exists(Sys.getenv("AVESDATA_PATH"))) {
  clootl::get_avesdata_repo(path = raw_path())
}
.check <- pr_get_tree(spp_data, source = "clootl", n_tree = 1)
if (length(.check$unmatched) > 0) {
  message("Removing ", length(.check$unmatched), " species absent from eBird taxonomy: ",
          paste(.check$unmatched, collapse = ", "))
  spp_data <- setdiff(spp_data, .check$unmatched)
}
stopifnot(length(.check$unmatched) <= 1)
got    <- pr_get_tree(spp_data, source = "clootl", n_tree = 100, cache = TRUE)
trees  <- got$tree
clootl_version <- pr_cite_tree(got, format = "text")

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
tree_samp <- sample(trees_pruned, 50)

make_A <- function(tree) {
  tree$tip.label <- gsub(" ", "_", tree$tip.label)
  if (!ape::is.ultrametric(tree)) tree <- phytools::force.ultrametric(tree, method = "nnls")
  inv <- inverseA(tree, nodes = "TIPS", scale = TRUE)
  A   <- solve(inv$Ainv); rownames(A) <- rownames(inv$Ainv); A
}

# --- log(body-mass) model (same structure as wing; mass-scale intercept prior) ---
priors_mass <- c(prior(normal(3, 1.5), class = Intercept),   # log grams ~ 3
                 prior(normal(0, 10),  class = b),
                 prior(cauchy(0, 1),   class = sd),
                 prior(cauchy(0, 1),   class = sigma))

fit_mass <- function(tree) {
  A <- make_A(tree)
  brm(ln_body_mass ~ 1 + Sex + scaled_yr + scaled_lat +
        (1 + scaled_yr || spp) + (1 | gr(species_name, cov = A)),
      data = dat, data2 = list(A = A),
      family = gaussian(), prior = priors_mass,
      iter = 1500, warmup = 1000, chains = 1, cores = 1,
      control = list(max_treedepth = 12, adapt_delta = 0.95),
      seed = 20240101)
}

plan(multisession)
mass_fits <- future_lapply(tree_samp, fit_mass, future.seed = TRUE)

# --- Rubin's rules pooling across trees (identical to atlantic_parallel.R) ----
pool_rubin <- function(fits, pars = c("b_Intercept","b_SexMale","b_scaled_yr","b_scaled_lat")) {
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
mass_rubin_summary <- pool_rubin(mass_fits)
print(mass_rubin_summary)

# phylogenetic signal (proportion of variance), pooled across the 50 fits
mass_sig <- sapply(mass_fits, function(f) {
  hypothesis(f,
    "sd_species_name__Intercept^2 / (sd_species_name__Intercept^2 + sd_spp__Intercept^2 + sd_spp__scaled_yr^2 + sigma^2) = 0",
    class = NULL)$hypothesis$Estimate
})
print(round(c(mean = mean(mass_sig), lwr = quantile(mass_sig, .025), upr = quantile(mass_sig, .975)), 3))

n_obs <- nobs(mass_fits[[1]]); n_spp_mod <- summary(mass_fits[[1]])$ngrps$species_name
save(mass_fits, mass_rubin_summary, mass_sig, clootl_version,
     file = out_path("models", "brm_mass_multiphylo.rda"))

# ============================================================================
# Model-free quartile comparison (mirrors the wing/bill descriptive check) +
# results summary, so the body-mass trend can be read directly against wing.
# ============================================================================
pq <- passer90 %>% filter(Year <= 2006 | Year >= 2013)
pq$period <- ifelse(pq$Year <= 2006, "1990-2006", "2013-2018")
mass.by.spp <- pq %>%
  group_by(spp = Binomial, period) %>%
  summarise(mean.lnmass = mean(ln_body_mass, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(trend = mean.lnmass - lag(mean.lnmass)) %>% ungroup()
q_dec <- sum(mass.by.spp$trend < 0, na.rm = TRUE)
q_inc <- sum(mass.by.spp$trend > 0, na.rm = TRUE)

pe <- function(par, col) round(mass_rubin_summary[[col]][mass_rubin_summary$par == par], 3)
results <- list(
  generated = as.character(Sys.time()),
  n_records = nrow(passer90), n_mass_records = n_obs,
  n_species = length(unique(passer90$Binomial)), n_species_model = as.integer(n_spp_mod),
  year_range = range(passer90$Year),
  rubin = mass_rubin_summary,
  phylo_signal = c(mean = mean(mass_sig), lwr = quantile(mass_sig, .025), upr = quantile(mass_sig, .975)),
  quartile = c(decreased = q_dec, increased = q_inc),
  year_beta = pe("b_scaled_yr", "estimate"),
  year_ci   = c(pe("b_scaled_yr", "lower"), pe("b_scaled_yr", "upper"))
)
saveRDS(results, out_path("mass_results.rds"))

md <- c(
  "# log(body mass) trend — 50-tree Rubin-pooled model (mirrors atlantic_parallel.R)",
  paste0("_generated ", results$generated, " · ", results$n_mass_records,
         " body-mass records across ", results$n_species_model, " species (sample: ",
         results$n_records, " records / ", results$n_species, " species), years ",
         results$year_range[1], "–", results$year_range[2], "_"),
  "", "## Pooled fixed effects (Wald–Rubin 95% CI)", "```",
  capture.output(print(mass_rubin_summary)), "```",
  "", paste0("Phylogenetic signal (proportion of among-species variance): ",
             sprintf("%.2f [%.2f, %.2f]", results$phylo_signal[1],
                     results$phylo_signal[2], results$phylo_signal[3])),
  "", paste0("Model-free quartile comparison (late vs early period): ",
             q_dec, " species decreased, ", q_inc, " increased in mean log mass."),
  "", paste0("**Year effect on log(body mass): β = ", results$year_beta,
             " [", results$year_ci[1], ", ", results$year_ci[2], "] per SD-year.**")
)
writeLines(md, out_path("mass_results.md"))
cat("\nWrote brm_mass_multiphylo.rda, mass_results.rds, mass_results.md\n")
