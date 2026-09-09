# _phylo_engine.R
# ---------------------------------------------------------------------------
# Shared phylogenetic-model engine for the revision (source() this from scripts).
#
# WHAT: fits phylogenetic mixed models with glmmTMB's `propto` covariance structure
#       (Williams, McGillycuddy, Drobniak, Bolker, Warton & Nakagawa 2025, bioRxiv
#       10.64898/2025.12.20.695312) across a sample of trees and pools with Rubin's
#       rules (Nakagawa & de Villemereuil 2019), replacing the brms 50-tree runs.
# WHY:  glmmTMB reproduces the published brms wing model to three decimals
#       (glmmtmb_validation_wing.R; GLMMTMB_ENGINE.md) at ~0.4 s per fit instead of
#       hours, so every model in the revision can carry phylogeny across 50 trees.
# TREES: phylo_A_list() returns the per-tree species CORRELATION matrices. By default
#       it uses the 50 matrices stored in the published brms fits (brm0_multiphylo.rda,
#       cached once to data/derived/phylo_A_50trees.rds) so all revision models use
#       the IDENTICAL trees as the published analysis. If a species set is not covered
#       it falls back to the clootl / prepR4pcm retrieval copied from
#       atlantic_parallel.R (same seed, same 50-of-100 sample).
# REQUIRES: glmmTMB >= 1.1.14 (CRAN; `propto` is a valid covstruct there).
# ---------------------------------------------------------------------------
suppressMessages({ library(glmmTMB); library(dplyr) })

# ABT (2018) -> current eBird/Clements names (identical map to atlantic_parallel.R)
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

#' Map ABT binomials (underscored) to the tree tip names (underscored eBird names).
phylo_species_name <- function(binomial) {
  nm <- gsub("_", " ", binomial)
  hit <- nm %in% names(ebird_synonyms); nm[hit] <- ebird_synonyms[nm[hit]]
  gsub(" ", "_", nm)
}

#' Build the A cache from the published brms fits (run once; ~5 s).
phylo_A_cache_from_brms <- function(rda = out_path("models", "brm0_multiphylo.rda"),
                                    cache = derived_path("phylo_A_50trees.rds")) {
  e <- new.env(); suppressWarnings(load(rda, envir = e))
  A_list <- lapply(e$fits, function(f) { A <- f$data2$A; A[order(rownames(A)), order(colnames(A))] })
  stopifnot(length(A_list) == 50, all(abs(diag(A_list[[1]]) - 1) < 1e-6))
  obj <- list(generated = Sys.time(), source = basename(rda), n_trees = length(A_list),
              species = rownames(A_list[[1]]), A = A_list,
              note = "Species correlation matrices (inverseA scale=TRUE, force.ultrametric nnls) of the 50 clootl trees sampled in atlantic_parallel.R (seed 20240101).")
  saveRDS(obj, cache); obj
}

#' Retrieve per-tree correlation matrices for a species set (tip names, underscored).
phylo_A_list <- function(species, n_trees = 50, cache = derived_path("phylo_A_50trees.rds"), seed = 20240101) {
  if (!file.exists(cache) && file.exists(out_path("models", "brm0_multiphylo.rda"))) phylo_A_cache_from_brms(cache = cache)
  if (file.exists(cache)) {
    obj <- readRDS(cache)
    missing <- setdiff(species, obj$species)
    if (!length(missing)) {
      return(lapply(obj$A[seq_len(min(n_trees, obj$n_trees))], function(A) A[species, species, drop = FALSE]))
    }
    message("[phylo] ", length(missing), " species not in the cached trees (", paste(head(missing, 5), collapse = ", "),
            "); rebuilding from clootl.")
  }
  # ---- fallback: identical retrieval to atlantic_parallel.R ----------------------
  suppressMessages({ library(ape); library(prepR4pcm); library(phytools); library(MCMCglmm) })
  set.seed(seed)
  aves <- raw_path("AvesDataLite-main")
  if (!nzchar(Sys.getenv("AVESDATA_PATH")) || !dir.exists(Sys.getenv("AVESDATA_PATH"))) {
    if (dir.exists(aves)) Sys.setenv(AVESDATA_PATH = aves) else clootl::get_avesdata_repo(path = raw_path())
  }
  spp_space <- gsub("_", " ", species)
  chk <- pr_get_tree(spp_space, source = "clootl", n_tree = 1)
  if (length(chk$unmatched)) { message("[phylo] dropping unmatched: ", paste(chk$unmatched, collapse = ", ")); spp_space <- setdiff(spp_space, chk$unmatched) }
  got <- pr_get_tree(spp_space, source = "clootl", n_tree = 100, cache = TRUE)
  keep <- gsub(" ", "_", spp_space)
  trees <- lapply(got$tree, function(t) { t$tip.label <- gsub(" ", "_", t$tip.label); ape::keep.tip(t, intersect(keep, t$tip.label)) })
  class(trees) <- "multiPhylo"
  samp <- sample(trees, n_trees)
  lapply(samp, function(tree) {
    if (!ape::is.ultrametric(tree)) tree <- phytools::force.ultrametric(tree, method = "nnls")
    inv <- inverseA(tree, nodes = "TIPS", scale = TRUE); A <- solve(inv$Ainv); rownames(A) <- colnames(A) <- rownames(inv$Ainv)
    A <- as.matrix(A); A[order(rownames(A)), order(colnames(A))]
  })
}

#' Fit one phylogenetic glmmTMB model. `formula` is the non-phylogenetic part; the
#' phylogenetic species intercept propto(0 + <species_col> | g, A) is appended.
fit_phylo_glmmtmb <- function(formula, data, A, species_col = "species_name", dispformula = ~1,
                              family = gaussian(), REML = TRUE, control = glmmTMBControl(), ...) {
  stopifnot(species_col %in% names(data))
  sp <- as.character(data[[species_col]])
  miss <- setdiff(unique(sp), rownames(A)); if (length(miss)) stop("species not in A: ", paste(head(miss), collapse = ", "))
  data <- data[sp %in% rownames(A), , drop = FALSE]
  data[[species_col]] <- factor(as.character(data[[species_col]]), levels = rownames(A))
  data$g <- factor(1)
  f <- update(formula, as.formula(paste(". ~ . + propto(0 +", species_col, "| g, A)")))
  environment(f) <- list2env(list(A = A), parent = environment(formula))
  glmmTMB(f, data = data, dispformula = dispformula, family = family, REML = REML, control = control, ...)
}

#' Tidy one fit: fixed effects (cond + disp), variance components, phylo proportion.
tidy_phylo_fit <- function(fit, tree = NA_integer_) {
  co <- summary(fit)$coefficients
  fe <- do.call(rbind, lapply(names(co), function(comp) if (!is.null(co[[comp]]) && nrow(co[[comp]]))
    data.frame(component = comp, par = rownames(co[[comp]]), estimate = co[[comp]][, 1], se = co[[comp]][, 2], row.names = NULL)))
  vc <- VarCorr(fit)$cond
  sds <- unlist(lapply(names(vc), function(nm) { s <- attr(vc[[nm]], "stddev"); setNames(s, paste0(nm, ":", names(s))) }))
  sd_phylo <- unname(sds[grep("^g:", names(sds))][1])
  sig <- tryCatch(sigma(fit), error = function(e) NA_real_)
  other <- sds[!grepl("^g:", names(sds))]
  list(fixed = cbind(tree = tree, fe),
       varcomp = data.frame(tree = tree, sd_phylo = sd_phylo, sigma = sig,
                            phylo_prop = sd_phylo^2 / (sd_phylo^2 + sum(other^2) + ifelse(is.na(sig), 0, sig^2)),
                            t(as.data.frame(as.list(other), check.names = FALSE)) |> (\(m) as.data.frame(t(m)))(),
                            converged = isTRUE(fit$sdr$pdHess), row.names = NULL))
}

#' Rubin's rules over per-tree fixed-effect tables.
pool_rubin_df <- function(fixed_tbl) {
  fixed_tbl %>% group_by(component, par) %>%
    summarise(m = n(), qbar = mean(estimate), ubar = mean(se^2), b = ifelse(n() > 1, var(estimate), 0), .groups = "drop") %>%
    mutate(se = sqrt(ubar + (1 + 1/m) * b), estimate = qbar, lower = qbar - 1.96 * se, upper = qbar + 1.96 * se,
           z = qbar / se) %>% select(component, par, m, estimate, se, lower, upper, z, ubar, between_tree_var = b)
}

#' Fit across trees and pool. Returns list(pooled, per_tree_fixed, varcomp, fits (optional), secs).
run_phylo_trees <- function(formula, data, A_list, keep_fits = FALSE, verbose = TRUE, ...) {
  t0 <- Sys.time(); fits <- vector("list", length(A_list)); tidy <- vector("list", length(A_list))
  for (i in seq_along(A_list)) {
    fit <- fit_phylo_glmmtmb(formula, data, A_list[[i]], ...)
    tidy[[i]] <- tidy_phylo_fit(fit, tree = i); if (keep_fits) fits[[i]] <- fit
    if (verbose && (i %% 10 == 0 || i == length(A_list))) message("[phylo] tree ", i, "/", length(A_list), " ", round(as.numeric(Sys.time() - t0, units = "secs")), " s")
  }
  fixed <- do.call(rbind, lapply(tidy, `[[`, "fixed")); varcomp <- do.call(rbind, lapply(tidy, `[[`, "varcomp"))
  list(pooled = pool_rubin_df(fixed), per_tree_fixed = fixed, varcomp = varcomp,
       varcomp_summary = varcomp %>% summarise(across(where(is.numeric) & !tree, list(mean = mean, lo = ~quantile(.x, .025), hi = ~quantile(.x, .975)))),
       fits = if (keep_fits) fits else NULL, n_trees = length(A_list), secs = as.numeric(Sys.time() - t0, units = "secs"),
       glmmTMB_version = as.character(packageVersion("glmmTMB")))
}

#' Species-specific slopes for a covariate from one fit: fixed + BLUP of the (spp) slope.
species_slopes_phylo <- function(fit, slope_term = "scaled_yr", spp_col = "spp") {
  re <- ranef(fit, condVar = TRUE)$cond[[spp_col]]
  if (is.null(re) || !slope_term %in% names(re)) stop("no random slope '", slope_term, "' for ", spp_col)
  cv <- attr(re, "condVar"); j <- match(slope_term, names(re))
  se_blup <- if (!is.null(cv)) sqrt(cv[j, j, ]) else NA_real_
  fe <- fixef(fit)$cond[[slope_term]]; se_fe <- summary(fit)$coefficients$cond[slope_term, 2]
  data.frame(spp = rownames(re), ranef = re[[slope_term]], se_ranef = se_blup, slope = fe + re[[slope_term]],
             se_total = sqrt(se_fe^2 + se_blup^2), row.names = NULL) %>%
    mutate(lo = slope - 1.96 * se_total, hi = slope + 1.96 * se_total)
}
