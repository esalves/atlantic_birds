# atlantic_phylo_slopes.R
# ---------------------------------------------------------------------------
# WHAT: extends the glmmTMB phylogenetic engine (_phylo_engine.R) with a SECOND
# propto() covariance block so the species year-SLOPE, not just the species
# INTERCEPT, carries a Brownian (phylogenetic) structure -- Gustavo's review
# comment that the response to change can itself be an extended phenotype.
#
# Design columns of the slope block are species indicators x scaled_yr; its A
# matrix (Asl) is a copy of the tree correlation matrix A whose dimnames are
# relabelled to those design-column names, so:
#   ... + propto(0 + species_name | g, A) + propto(0 + species_name:scaled_yr | g, Asl)
# VarCorr(fit)$cond names these "g" (phylo intercept) and "g.1" (phylo slope);
# ranef(fit)$cond$g comes back as one 1 x 2S row (species intercepts, then
# species:scaled_yr slopes, confirmed empirically -- see tidy_pslope_fit()).
# _phylo_engine.R's tidy_phylo_fit()/fit_phylo_glmmtmb() are NOT extended (other
# scripts source that file concurrently); this script writes its own fitting
# and tidying functions instead.
#
# MODELS (P-ladder), all with a phylogenetic species INTERCEPT (as published):
#   P0  current published structure: (1 + scaled_yr || spp) [iid int + slope]
#       + propto(0 + species_name | g, A)                    [phylo intercept]
#   P1  P0 + propto(0 + species_name:scaled_yr | g, Asl)      [+ phylo slope]
#       (iid slope kept: partitions slope variance into phylo + iid)
#   P2  (1 | spp) [iid intercept only, iid slope DROPPED]
#       + propto(0 + species_name | g, A) + propto(0 + species_name:scaled_yr | g, Asl)
#       (phylogenetic slope only)
#
# RESPONSES: wing M0, wing M3 (REQUIRED; M3 reproduction checked against
#   output/controlled_wing_phylo_results.rds to 3 dp before proceeding), plus
#   two data are readily available for -- the E4b isometry contrast
#   (iso = lnmass - 3*lnwing, referee_reruns.R's M3-adjustment spec) and ln
#   body mass under the wing M3 spec, both on referee_reruns.R's `shared`
#   sample (wing + mass + climate complete, 7,577 records after M3 guards).
#
# COST: developed/smoke-tested on 2 trees; the 50-tree job is run ONCE from
# the repo root, logged to the session scratchpad (see REPORT.md for the log
# path). Re-runs of this script read the saved .rds instead of refitting.
#
# OUTPUT: output/phylo_slopes/phylo_slopes_results.rds, REPORT.md (written by
#   a short follow-up block that reads the .rds -- see the bottom of this file).
# RUN:  Rscript Analysis/scripts/atlantic_phylo_slopes.R --trees 50
#       Rscript Analysis/scripts/atlantic_phylo_slopes.R --trees 2   (dev/smoke)
# ---------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(glmmTMB) })
options(width = 140)
set.seed(20240101)

args <- commandArgs(trailingOnly = TRUE)
get_flag <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
N_TREES <- as.integer(get_flag("--trees", "2"))
message(sprintf("[phylo_slopes] N_TREES = %d", N_TREES))

# --- path helpers (same convention as the other scripts) -------------------
.find_analysis_dir <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:10) {
    if (dir.exists(file.path(d, "data", "derived")) && dir.exists(file.path(d, "scripts"))) return(d)
    if (dir.exists(file.path(d, "Analysis", "data", "derived"))) return(normalizePath(file.path(d, "Analysis"), winslash = "/"))
    parent <- dirname(d); if (identical(parent, d)) break; d <- parent
  }
  stop("Could not locate Analysis/. Run from inside the atlantic_birds repo.")
}
ANALYSIS_DIR <- .find_analysis_dir()
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
script_path  <- function(...) file.path(ANALYSIS_DIR, "scripts", ...)
OUT_DIR <- out_path("phylo_slopes"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source(script_path("_phylo_engine.R"))   # phylo_species_name(), phylo_A_list()

# ===========================================================================
# 1. Data: wing frame (trimmed copy of atlantic_parallel_controlled.R section 1)
# ===========================================================================
load(derived_path("passer90.rda"))
load(derived_path("passer90_allsex.rda"))
scaling <- attr(passer90, "scaling")
SD_YEAR <- unname(scaling$year["scale"]); stopifnot(abs(SD_YEAR - 5.02) < 0.01)
ALT_CENTRE <- unname(scaling$alt["center"]); ALT_SCALE <- unname(scaling$alt["scale"])

haversine_km <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180; a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * 6371 * asin(sqrt(pmin(1, a)))
}
alt_sites <- passer90_allsex %>% filter(!is.na(Altitude), !is.na(Latitude_decimal_degrees)) %>%
  group_by(Latitude_decimal_degrees, Longitude_decimal_degrees) %>% summarise(alt_site = median(Altitude), .groups = "drop")
na_sites <- passer90_allsex %>% filter(is.na(Altitude), !is.na(Latitude_decimal_degrees)) %>%
  distinct(Latitude_decimal_degrees, Longitude_decimal_degrees, Municipality, Locality, Main_researcher)
na_sites$alt_imp <- NA_real_; na_sites$dist_km <- NA_real_
for (i in seq_len(nrow(na_sites))) {
  dk <- haversine_km(na_sites$Latitude_decimal_degrees[i], na_sites$Longitude_decimal_degrees[i],
                     alt_sites$Latitude_decimal_degrees, alt_sites$Longitude_decimal_degrees)
  k <- which.min(dk); na_sites$alt_imp[i] <- alt_sites$alt_site[k]; na_sites$dist_km[i] <- dk[k]
}
alt_lookup <- na_sites %>% select(Latitude_decimal_degrees, Longitude_decimal_degrees, alt_imp, dist_km) %>% distinct()

prep_frame <- function(d) {
  d %>% left_join(alt_lookup, by = c("Latitude_decimal_degrees", "Longitude_decimal_degrees")) %>%
    mutate(scaled_yr = as.numeric(scaled_yr), scaled_lat = as.numeric(scaled_lat), scaled_lon = as.numeric(scaled_lon),
           alt_imputed = is.na(Altitude), alt_filled = ifelse(alt_imputed, alt_imp, Altitude),
           scaled_alt = ifelse(is.na(alt_filled), 0, (alt_filled - ALT_CENTRE) / ALT_SCALE),
           Sex = factor(Sex, levels = c("Female", "Male")), spp0 = factor(Binomial),
           src = factor(Main_researcher), site = factor(Municipality),
           ind = factor(ifelse(is.na(Ring), paste0("unringed_", ID_ABT), paste(Ring, Binomial, sep = "|"))),
           wing_col = factor(wing_col, levels = c("right", "left", "generic")),
           season = factor(ifelse(is.na(season), "unrecorded", as.character(season)), levels = c("DJF", "MAM", "JJA", "SON", "unrecorded")),
           molt3 = factor(ifelse(is.na(Molt), "unrecorded", Molt), levels = c("no", "Yes", "unrecorded")),
           has_wing = !is.na(conc.wing.length))
}
known <- prep_frame(passer90)
w_all <- known %>% filter(has_wing) %>% droplevels()
stopifnot(nrow(w_all) == 8478)

# --- relabel to tree tip names (species_name / spp), drop species absent from the trees ----
cache_file <- derived_path("phylo_A_50trees.rds")
tree_species <- readRDS(cache_file)$species
attach_phylo_names <- function(d) {
  d$species_name <- phylo_species_name(d$Binomial)
  drop <- setdiff(unique(d$species_name), tree_species)
  if (length(drop)) message("[phylo_slopes] dropping ", length(drop), " species absent from the trees: ", paste(drop, collapse = ", "))
  d %>% filter(!species_name %in% drop) %>% mutate(spp = factor(species_name)) %>% droplevels()
}
w_all_p <- attach_phylo_names(w_all)
stopifnot(nrow(w_all_p) == nrow(w_all))   # no species lost on the wing sample (matches controlled tier's message)

# ===========================================================================
# 2. Generic P0/P1/P2 fitting + tidying engine
# ===========================================================================
#' Fit one P-variant of one model on one tree.
#' iid_re: "(1 + scaled_yr || spp)" [P0/P1] or "(1 | spp)" [P2]
#' phylo_slope: FALSE [P0] or TRUE [P1/P2]
fit_pslope <- function(response, fixed, iid_re, extra_re, data, A, phylo_slope, REML = TRUE, control = glmmTMBControl()) {
  sp <- as.character(data$species_name)
  stopifnot(all(sp %in% rownames(A)))
  data <- data[sp %in% rownames(A), , drop = FALSE]
  data$species_name <- factor(as.character(data$species_name), levels = rownames(A))
  data$g <- factor(1)
  re_parts <- c(iid_re, if (nzchar(extra_re)) extra_re, "propto(0 + species_name | g, A)")
  if (phylo_slope) re_parts <- c(re_parts, "propto(0 + species_name:scaled_yr | g, Asl)")
  f <- as.formula(paste(response, "~ 1 +", fixed, "+", paste(re_parts, collapse = " + ")))
  Asl <- NULL
  if (phylo_slope) { Asl <- A; dimnames(Asl) <- rep(list(paste0("species_name", rownames(A), ":scaled_yr")), 2) }
  environment(f) <- list2env(list(A = A, Asl = Asl), parent = globalenv())
  glmmTMB(f, data = data, REML = REML, control = control)
}

#' Tidy one fit: fixed effects + variance components (own naming; does NOT use
#' _phylo_engine.R's tidy_phylo_fit(), which would lump g.1 into "other").
tidy_pslope_fit <- function(fit, tree) {
  co <- summary(fit)$coefficients$cond
  fe <- data.frame(tree = tree, par = rownames(co), estimate = co[, 1], se = co[, 2], row.names = NULL)
  vc <- VarCorr(fit)$cond
  # "g"/"g.1" (propto blocks) carry ONE shared variance broadcast across every species
  # (dimnames are per-species design columns, not "(Intercept)"/"scaled_yr"): take element 1.
  phylo_sd <- function(grp) if (grp %in% names(vc)) unname(attr(vc[[grp]], "stddev")[1]) else NA_real_
  # "spp" (iid block) has genuinely different intercept/slope SDs: look up by name.
  iid_sd <- function(nm) { if (!"spp" %in% names(vc)) return(NA_real_); s <- attr(vc[["spp"]], "stddev"); if (nm %in% names(s)) unname(s[[nm]]) else NA_real_ }
  other_grps <- setdiff(names(vc), c("g", "g.1", "spp"))
  other_var <- if (length(other_grps)) sum(vapply(other_grps, function(g) unname(attr(vc[[g]], "stddev")[1])^2, numeric(1))) else 0
  data.frame(tree = tree,
             sd_phylo_int = phylo_sd("g"), sd_phylo_slope = phylo_sd("g.1"),
             sd_iid_int = iid_sd("(Intercept)"), sd_iid_slope = iid_sd("scaled_yr"),
             sigma = tryCatch(sigma(fit), error = function(e) NA_real_),
             other_grps = paste(other_grps, collapse = ","), other_var = other_var,
             logLik = as.numeric(logLik(fit)), pdHess = isTRUE(fit$sdr$pdHess), row.names = NULL)
}

#' Species slope decomposition under P1/P2: slope_i = beta_yr + u_iid_i + u_phylo_i.
species_slopes_pslope <- function(fit, A, has_iid_slope) {
  fe <- fixef(fit)$cond[["scaled_yr"]]
  spp <- rownames(A)
  re_spp <- tryCatch(ranef(fit)$cond$spp, error = function(e) NULL)
  u_iid <- setNames(rep(0, length(spp)), spp)
  if (has_iid_slope && !is.null(re_spp) && "scaled_yr" %in% colnames(re_spp)) {
    u <- setNames(re_spp[["scaled_yr"]], rownames(re_spp)); u_iid[names(u)] <- u
  }
  re_g <- ranef(fit)$cond$g; cn <- colnames(re_g)
  slope_cols <- grep(":scaled_yr$", cn)
  nm <- sub("^species_name(.*):scaled_yr$", "\\1", cn[slope_cols])
  u_phylo <- setNames(as.numeric(re_g[1, slope_cols]), nm)[spp]
  data.frame(spp = spp, u_iid = u_iid[spp], u_phylo = u_phylo, fixed = fe,
             slope = fe + u_iid[spp] + u_phylo, row.names = NULL)
}

#' Run P0/P1/P2 for one model spec across A_list, pool with Rubin's rules, return everything.
run_pladder <- function(model_name, response, fixed, iid_re_p0, extra_re, data, A_list, sd_year, mean_resp = NA, is_log = FALSE) {
  n <- length(A_list)
  fixed_rows <- list(P0 = list(), P1 = list(), P2 = list())
  vc_rows    <- list(P0 = list(), P1 = list(), P2 = list())
  slopes_P1  <- vector("list", n)
  fits_keep  <- list(P0 = vector("list", n), P1 = vector("list", n))   # for LRT (logLik already in vc_rows; keep species slopes only)
  t0 <- Sys.time()
  for (i in seq_len(n)) {
    A <- A_list[[i]]
    fit0 <- fit_pslope(response, fixed, iid_re_p0, extra_re, data, A, phylo_slope = FALSE)
    fit1 <- fit_pslope(response, fixed, iid_re_p0, extra_re, data, A, phylo_slope = TRUE)
    fit2 <- fit_pslope(response, fixed, "(1 | spp)", extra_re, data, A, phylo_slope = TRUE)
    tf0 <- tidy_pslope_fit(fit0, i); tf1 <- tidy_pslope_fit(fit1, i); tf2 <- tidy_pslope_fit(fit2, i)
    fixed_rows$P0[[i]] <- data.frame(tree = i, par = rownames(summary(fit0)$coefficients$cond),
                                     estimate = summary(fit0)$coefficients$cond[, 1], se = summary(fit0)$coefficients$cond[, 2], row.names = NULL)
    fixed_rows$P1[[i]] <- data.frame(tree = i, par = rownames(summary(fit1)$coefficients$cond),
                                     estimate = summary(fit1)$coefficients$cond[, 1], se = summary(fit1)$coefficients$cond[, 2], row.names = NULL)
    fixed_rows$P2[[i]] <- data.frame(tree = i, par = rownames(summary(fit2)$coefficients$cond),
                                     estimate = summary(fit2)$coefficients$cond[, 1], se = summary(fit2)$coefficients$cond[, 2], row.names = NULL)
    vc_rows$P0[[i]] <- tf0; vc_rows$P1[[i]] <- tf1; vc_rows$P2[[i]] <- tf2
    slopes_P1[[i]] <- cbind(tree = i, species_slopes_pslope(fit1, A, has_iid_slope = TRUE))
    if (i %% 10 == 0 || i == n) message(sprintf("[phylo_slopes] %s tree %d/%d (%.0fs)", model_name, i, n, as.numeric(Sys.time() - t0, units = "secs")))
  }
  pool <- function(rows) {
    ft <- bind_rows(rows); ft$component <- "cond"
    p <- pool_rubin_df(ft, on_drop = "warn")
    if (is_log) {
      p$per_decade <- 100 * (exp(p$estimate * 10 / sd_year) - 1)       # % / decade (log-scale response)
      p$pct_per_decade <- p$per_decade
    } else {
      p$per_decade <- p$estimate * 10 / sd_year                        # mm / decade
      p$pct_per_decade <- if (!is.na(mean_resp)) p$per_decade / mean_resp * 100 else NA_real_
    }
    p
  }
  pooled <- list(P0 = pool(fixed_rows$P0), P1 = pool(fixed_rows$P1), P2 = pool(fixed_rows$P2))
  vc <- list(P0 = bind_rows(vc_rows$P0), P1 = bind_rows(vc_rows$P1), P2 = bind_rows(vc_rows$P2))
  H2_slope <- with(vc$P1, sd_phylo_slope^2 / (sd_phylo_slope^2 + sd_iid_slope^2))
  lrt <- data.frame(tree = seq_len(n), logLik_P0 = vc$P0$logLik, logLik_P1 = vc$P1$logLik) %>%
    mutate(stat = pmax(0, 2 * (logLik_P1 - logLik_P0)), p_boundary = 0.5 * pchisq(stat, df = 1, lower.tail = FALSE))
  list(model = model_name, n_trees = n, fixed = fixed_rows, pooled = pooled, varcomp = vc,
       H2_slope = H2_slope, H2_slope_summary = c(mean = mean(H2_slope, na.rm = TRUE), lo = unname(quantile(H2_slope, .025, na.rm = TRUE)), hi = unname(quantile(H2_slope, .975, na.rm = TRUE))),
       lrt = lrt, lrt_summary = c(mean_stat = mean(lrt$stat), median_p = median(lrt$p_boundary), pct_p_lt_05 = mean(lrt$p_boundary < 0.05)),
       converged = c(P0 = sum(vc$P0$pdHess), P1 = sum(vc$P1$pdHess), P2 = sum(vc$P2$pdHess)),
       species_slopes_P1 = bind_rows(slopes_P1), sd_year = sd_year, secs = as.numeric(Sys.time() - t0, units = "secs"))
}

FX_M0 <- "Sex + scaled_yr + scaled_lat"
FX_M3 <- "Sex + scaled_yr + scaled_lat + wing_col + scaled_lon + scaled_alt + season + molt3"
RE_M3_EXTRA <- "(1 | src) + (1 | site) + (1 | ind)"
RE_P0 <- "(1 + scaled_yr || spp)"

A_wing <- phylo_A_list(sort(unique(as.character(w_all_p$spp))), n_trees = N_TREES)
mean_wing <- mean(w_all_p$conc.wing.length)

message("[phylo_slopes] === wing M0 ===")
res_wingM0 <- run_pladder("wing_M0", "conc.wing.length", FX_M0, RE_P0, "", w_all_p, A_wing, SD_YEAR, mean_wing, is_log = FALSE)
message("[phylo_slopes] === wing M3 ===")
res_wingM3 <- run_pladder("wing_M3", "conc.wing.length", FX_M3, RE_P0, RE_M3_EXTRA, w_all_p, A_wing, SD_YEAR, mean_wing, is_log = FALSE)

# --- REPRODUCTION CHECK: wing M0/M3 under P0 must match controlled_wing_phylo_results.rds ----
repro <- NULL
val_file <- out_path("controlled_wing_phylo_results.rds")
if (file.exists(val_file)) {
  ref <- readRDS(val_file)
  get_ref <- function(m) ref$before_after %>% filter(model == m, term == "scaled_yr") %>% select(estimate, ci_lo, ci_hi)
  get_new <- function(res) res$pooled$P0 %>% filter(par == "scaled_yr") %>% select(estimate, lower, upper)
  r0 <- get_ref("M0"); n0 <- get_new(res_wingM0)
  r3 <- get_ref("M3"); n3 <- get_new(res_wingM3)
  repro <- list(M0 = list(reference = r0, new_P0 = n0, match_3dp = all(abs(unlist(r0) - unlist(n0)) < 5e-4)),
               M3 = list(reference = r3, new_P0 = n3, match_3dp = all(abs(unlist(r3) - unlist(n3)) < 5e-4)))
  message(sprintf("[repro] M0: ref=%.4f new=%.4f | M3: ref=%.4f new=%.4f | match_3dp: M0=%s M3=%s",
                  r0$estimate, n0$estimate, r3$estimate, n3$estimate, repro$M0$match_3dp, repro$M3$match_3dp))
} else message("[repro] controlled_wing_phylo_results.rds not found; skipping reproduction check")

# ===========================================================================
# 3. Isometry (E4b) and ln body mass (M3 spec), referee_reruns.R `shared` sample
# ===========================================================================
res_iso <- NULL; res_mass <- NULL
pc <- tryCatch(readRDS(derived_path("passer90_climate.rds")), error = function(e) NULL)
if (!is.null(pc)) {
  shared0 <- pc %>% mutate(
    scaled_yr = as.numeric(scaled_yr), scaled_lat = as.numeric(scaled_lat), scaled_lon = as.numeric(scaled_lon), scaled_alt = as.numeric(scaled_alt),
    Sex = factor(Sex, levels = c("Female", "Male")), src = factor(Main_researcher), site = factor(Municipality),
    season = factor(season, levels = c("DJF", "MAM", "JJA", "SON")), molt = factor(ifelse(is.na(Molt), "unknown", as.character(Molt))),
    ind = factor(paste(Ring, Binomial, sep = "_")), wing = conc.wing.length, lnwing = log(conc.wing.length), lnmass = ln_body_mass,
    wing_col = factor(wing_col, levels = c("right", "left", "generic")),
    molt3 = factor(ifelse(is.na(Molt), "unrecorded", Molt), levels = c("no", "Yes", "unrecorded")),
    species_name = phylo_species_name(Binomial)) %>%
    filter(!is.na(wing), !is.na(lnmass), !is.na(scaled_lon), !is.na(scaled_alt), !is.na(season)) %>%
    mutate(iso = lnmass - 3 * lnwing) %>% as.data.frame()
  drop_sp <- setdiff(unique(shared0$species_name), tree_species)
  shared <- shared0 %>% filter(!species_name %in% drop_sp) %>% mutate(spp = factor(species_name)) %>% droplevels()
  SD_YR_SHARED <- as.numeric(attr(pc$scaled_yr, "scaled:scale"))
  A_shared <- phylo_A_list(sort(unique(as.character(shared$spp))), n_trees = N_TREES)
  message(sprintf("[phylo_slopes] shared sample (E4b/mass): %d records, %d species", nrow(shared), nlevels(shared$spp)))

  FX_E4B <- "Sex + scaled_yr + scaled_lat + scaled_lon + scaled_alt + season + molt"
  RE_E4B_EXTRA <- "(1 | src) + (1 | ind)"
  message("[phylo_slopes] === isometry (E4b) ===")
  res_iso <- run_pladder("iso_E4b", "iso", FX_E4B, RE_P0, RE_E4B_EXTRA, shared, A_shared, SD_YR_SHARED, is_log = TRUE)

  message("[phylo_slopes] === ln body mass (wing M3 spec) ===")
  res_mass <- run_pladder("lnmass_M3", "lnmass", FX_M3, RE_P0, RE_M3_EXTRA, shared, A_shared, SD_YR_SHARED, is_log = TRUE)
}

# ===========================================================================
# 4. Save
# ===========================================================================
out <- list(generated = Sys.time(), n_trees = N_TREES,
           glmmTMB_version = as.character(packageVersion("glmmTMB")),
           reproduction_check = repro,
           wing_M0 = res_wingM0, wing_M3 = res_wingM3, iso_E4b = res_iso, lnmass_M3 = res_mass,
           sd_year_wing = SD_YEAR, mean_wing_mm = mean_wing,
           note = paste("P0 = published structure (iid species intercept+slope + phylogenetic intercept);",
                        "P1 = P0 + phylogenetic slope (iid slope kept, partitions slope variance);",
                        "P2 = phylogenetic slope only (iid slope dropped, iid intercept kept).",
                        "H2_slope = sd_phylo_slope^2 / (sd_phylo_slope^2 + sd_iid_slope^2) under P1, per tree, then mean + 2.5/97.5% across trees.",
                        "LRT P0 vs P1: REML logLik (same fixed effects both sides), boundary-corrected p = 0.5*pchisq(stat,1,lower=FALSE).",
                        "per_decade = mm/decade (wing) or %/decade = 100*(exp(estimate*10/SD_year)-1) (iso, lnmass)."))
saveRDS(out, file.path(OUT_DIR, "phylo_slopes_results.rds"))
message("[write] ", file.path(OUT_DIR, "phylo_slopes_results.rds"))
