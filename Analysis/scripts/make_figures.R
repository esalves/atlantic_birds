# make_figures.R
# ---------------------------------------------------------------------------
# WHAT: Two jobs, selected by the positional argument(s):
#
#   (a) "export"  -> writes the CSV files that Manuscript/make_figures.py reads
#                    (Analysis/output/figure_data/*.csv). The Python figure
#                    script cannot read .rds / .rda files (pyreadr is not
#                    installed), so every table it needs is flattened here:
#                    the record-level analytical sample (passer90.rda), the
#                    Phase 0 audit tables (audit_sites.rds, audit_sources.rds),
#                    the effect-scale constants (effect_scale.rds), the
#                    Phase 1 lme4 fast-tier results (controlled_wing_results.rds,
#                    controlled_wing_species_slopes.rds) and, when present, the
#                    Phase 1 glmmTMB phylogenetic tier (controlled_wing_phylo_
#                    results.rds, controlled_wing_phylo_species_slopes.rds) which
#                    make_figures.py promotes over the lme4 tier as the primary
#                    fitted line / species-slope layer (2026-09 revision; see
#                    GLMMTMB_ENGINE.md). Missing inputs are skipped with a
#                    message, never fabricated; make_figures.py then skips the
#                    figures that depend on them or falls back a tier.
#   (b) "noarthro" / "arthro" (or no mode, = both) -> regenerate the drmSEM
#                    result figures from the saved Analysis/output/*.rds WITHOUT
#                    refitting any model (reads drmsem_results_<mode>.rds, which
#                    store paths_raw, paths_sdx and the effect tables as plain
#                    data frames) and redraws the four figures per mode into
#                    Analysis/figures/. Only needs ggplot2 (+ scales, grid).
#                    When no mode is given the export step (a) also runs.
#
# WHY (2026-09, revision Phase 6; REVISION_PLAN.md §3 "Phase 6"): the redesigned
# Fig. 1 (shared vs period-unique sites), the within-species-centred Fig. 2 with
# the controlled (M3) fit, the species-slope caterpillar and the provenance
# supplementary figures all need columns that passer90_export.csv (Binomial,
# Year, cwl, bill width, lon, lat) does not carry, plus the audit and Phase 1
# tables. Rather than teaching the Python script to read R serialisations, this
# script writes plain CSVs; make_figures.py calls it (Rscript ... export) before
# drawing, and falls back to the CSVs already on disk when Rscript is absent.
#
# INPUTS:  data/derived/passer90.rda
#          output/audit_sites.rds, output/audit_sources.rds, output/effect_scale.rds
#          output/controlled_wing_results.rds, output/controlled_wing_species_slopes.rds     (lme4 Tier 1)
#          output/controlled_wing_phylo_results.rds,                                         (glmmTMB phylo tier,
#            output/controlled_wing_phylo_species_slopes.rds                                  if present; optional)
#          output/controlled_wing_brms_results.rds                                           (brms Tier 2-3, optional)
#          output/drmsem_results_<mode>.rds                      (drmSEM figures)
# OUTPUTS: output/figure_data/fig_records.csv                  record-level sample
#          output/figure_data/fig_map_sites.csv                coordinate sites x period
#          output/figure_data/fig_contributors_per_year.csv
#          output/figure_data/fig_records_per_contributor_year.csv
#          output/figure_data/fig_contributor_table.csv
#          output/figure_data/fig_wingcol_by_year.csv
#          output/figure_data/fig_controlled_before_after.csv        Phase 1 year terms (lme4)
#          output/figure_data/fig_species_slopes.csv                 Phase 1 species slopes (lme4; M3, M3_cc, M0)
#          output/figure_data/fig_controlled_before_after_phylo.csv  Phase 1 year terms (glmmTMB phylo, if present)
#          output/figure_data/fig_species_slopes_phylo.csv           Phase 1 species slopes (glmmTMB phylo, if present)
#          output/figure_data/fig_scalars.csv                  key = value constants
#          figures/drmsem_*_<mode>.png (+ Manuscript/images/fig-drmsem.png)
#
# RUN:  Rscript Analysis/scripts/make_figures.R export     # CSVs for make_figures.py only
#       Rscript Analysis/scripts/make_figures.R            # export + both drmSEM modes
#       Rscript Analysis/scripts/make_figures.R arthro     # one drmSEM mode only
#       (works from the repo root, Analysis/, or Analysis/scripts/)
# Session: R 4.6.0, ggplot2 4.0.3, dplyr 1.2.1, data.table 1.17.x.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(ggplot2))
set.seed(20260909)   # nothing stochastic below; recorded for the convention

# --- robust path resolution (repo layout: Analysis/{scripts,data,output,figures}) ---
# Locate the Analysis/ directory regardless of getwd() (script may be run from the
# repo root, Analysis/, or Analysis/scripts/), then build paths into its subfolders.
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
ANALYSIS_DIR <- .find_analysis_dir()                 # = .../atlantic_birds/Analysis
raw_path     <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
derived_path <- function(...) file.path(ANALYSIS_DIR, "data", "derived", ...)
out_path     <- function(...) file.path(ANALYSIS_DIR, "output", ...)
fig_path     <- function(...) file.path(ANALYSIS_DIR, "figures", ...)
script_path  <- function(...) file.path(ANALYSIS_DIR, "scripts", ...)
.repo      <- dirname(ANALYSIS_DIR)
OUTPUT_DIR <- out_path()
FIG_DIR    <- fig_path()
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# ============================================================================
# (a) EXPORT: flatten the rds/rda inputs of make_figures.py to CSV
# ============================================================================
export_figure_data <- function() {
  suppressPackageStartupMessages(library(dplyr))
  fd <- out_path("figure_data"); dir.create(fd, showWarnings = FALSE, recursive = TRUE)
  wcsv <- function(df, name) {
    df <- as.data.frame(df)
    # drop list-columns (none expected) and flatten 1-column matrices from scale()
    for (v in names(df)) if (is.matrix(df[[v]])) df[[v]] <- as.numeric(df[[v]])
    data.table::fwrite(df, file.path(fd, name))
    message(sprintf("  wrote %-40s %6d rows x %2d cols", name, nrow(df), ncol(df)))
  }
  skip <- function(f) { message("  SKIP ", basename(f), " not found (the figures that need it will be skipped)"); FALSE }
  scal <- list()
  add_scalar <- function(key, value, source) scal[[length(scal) + 1]] <<- data.frame(key = key, value = as.character(value), source = source)

  # --- record-level analytical sample -----------------------------------------
  f <- derived_path("passer90.rda")
  if (file.exists(f)) {
    load(f)   # -> passer90 (known-sex, live, 73 spp, 12,571 records)
    sc <- attr(passer90, "scaling")
    rec <- passer90 %>%
      transmute(ID_ABT, Binomial, Family, Sex, Year, month, season = as.character(season),
                cwl = conc.wing.length, wing_col, Bill_width.mm., Body_mass.g.,
                Main_researcher, Municipality, Locality,
                Longitude_decimal_degrees, Latitude_decimal_degrees, Altitude, Molt, Ring)
    wcsv(rec, "fig_records.csv")
    add_scalar("n_records", nrow(rec), "passer90.rda")
    add_scalar("n_species", n_distinct(rec$Binomial), "passer90.rda")
    add_scalar("n_wing_records", sum(!is.na(rec$cwl)), "passer90.rda")
    add_scalar("year_centre_scaling", unname(sc$year["center"]), "attr(passer90, 'scaling')$year")
    add_scalar("year_sd_scaling", unname(sc$year["scale"]), "attr(passer90, 'scaling')$year")
  } else skip(f)

  # --- Phase 0 audit tables -----------------------------------------------------
  f <- out_path("audit_sites.rds")
  if (file.exists(f)) {
    a <- readRDS(f)
    wcsv(a$map_sites, "fig_map_sites.csv")
    for (k in c("locality_wing", "municipality_wing", "coordinate_site_wing", "species_municipality_wing"))
      for (kk in names(a[[k]])) add_scalar(paste(k, kk, sep = "."), a[[k]][[kk]], "audit_sites.rds")
    add_scalar("n_wing_quartile", a$n_wing_quartile, "audit_sites.rds")
  } else skip(f)
  f <- out_path("audit_sources.rds")
  if (file.exists(f)) {
    a <- readRDS(f)
    wcsv(a$contributors_per_year, "fig_contributors_per_year.csv")
    wcsv(a$records_per_contributor_year_wing, "fig_records_per_contributor_year.csv")
    wcsv(a$contributor_table, "fig_contributor_table.csv")
    wcsv(a$wingcol_by_year, "fig_wingcol_by_year.csv")
    add_scalar("n_wing_contributors", a$n_wing_contributors, "audit_sources.rds")
    add_scalar("n_span_both_wing", a$n_span_both_wing, "audit_sources.rds")
    add_scalar("wing_records_from_spanners", a$wing_records_from_spanners, "audit_sources.rds")
    add_scalar("prop_right_early", a$prop_right_early, "audit_sources.rds")
    add_scalar("prop_generic_late", a$prop_generic_late, "audit_sources.rds")
    add_scalar("spanning_contributors", paste(a$spanning_contributors, collapse = ";"), "audit_sources.rds")
  } else skip(f)
  f <- out_path("effect_scale.rds")
  if (file.exists(f)) {
    e <- readRDS(f)
    add_scalar("mean_wing_mm", e$mean_wing_mm, "effect_scale.rds")
    add_scalar("sd_year_scaling", e$sd_year_scaling, "effect_scale.rds")
    add_scalar("published_wing_est", e$published$wing[["est"]], "effect_scale.rds (brm0_multiphylo Rubin pooled)")
    add_scalar("published_wing_lo", e$published$wing[["lo"]], "effect_scale.rds")
    add_scalar("published_wing_hi", e$published$wing[["hi"]], "effect_scale.rds")
  } else skip(f)

  # --- Phase 1 fast-tier (lme4) results -------------------------------------------
  f <- out_path("controlled_wing_results.rds")
  if (file.exists(f)) {
    r <- readRDS(f)
    wcsv(r$before_after, "fig_controlled_before_after.csv")
    add_scalar("controlled_mode", r$mode, "controlled_wing_results.rds")
    add_scalar("controlled_generated", format(r$generated), "controlled_wing_results.rds")
    add_scalar("controlled_sd_year", r$sd_year, "controlled_wing_results.rds")
    add_scalar("controlled_mean_wing_mm", r$mean_wing_mm, "controlled_wing_results.rds")
    add_scalar("controlled_session", paste(names(r$session), unlist(r$session), sep = "=", collapse = "; "), "controlled_wing_results.rds")
    add_scalar("decision_gate_scenario", r$decision_gate$scenario, "controlled_wing_results.rds")
  } else skip(f)
  f <- out_path("controlled_wing_species_slopes.rds")
  if (file.exists(f)) {
    s <- readRDS(f)
    wcsv(bind_rows(s$M3, s$M3_cc, s$M0), "fig_species_slopes.csv")
    for (m in names(s$summary)) for (k in names(s$summary[[m]]))
      add_scalar(paste("slopes", m, k, sep = "."), s$summary[[m]][[k]], "controlled_wing_species_slopes.rds")
  } else skip(f)

  # --- Phase 1 phylogenetic tier (glmmTMB propto, Rubin-pooled over trees) --------
  # Default engine per REVISION_PLAN.md / GLMMTMB_ENGINE.md (Williams et al. 2025):
  # glmmTMB propto reproduces the published 50-tree brms wing model to 2-3 decimals
  # in seconds instead of hours (glmmtmb_validation_wing.R). When this file exists,
  # make_figures.py promotes it over the lme4 fast tier as the primary fitted line /
  # species-slope layer, keeping lme4 as a lighter non-phylogenetic comparison layer.
  # `$n_trees` in the file (not hard-coded here) is what actually appears in the
  # panel labels, so a partial/dev run is never mislabelled as the full 50-tree run.
  f <- out_path("controlled_wing_phylo_results.rds")
  if (file.exists(f)) {
    r <- readRDS(f)
    wcsv(r$before_after, "fig_controlled_before_after_phylo.csv")
    add_scalar("phylo_available", TRUE, "controlled_wing_phylo_results.rds")
    add_scalar("phylo_mode", r$mode, "controlled_wing_phylo_results.rds")
    add_scalar("phylo_engine_name", r$engine$name, "controlled_wing_phylo_results.rds")
    add_scalar("phylo_engine_version", r$engine$version, "controlled_wing_phylo_results.rds")
    add_scalar("phylo_n_trees", r$n_trees, "controlled_wing_phylo_results.rds")
    add_scalar("phylo_generated", format(r$generated), "controlled_wing_phylo_results.rds")
    add_scalar("phylo_sd_year", r$sd_year, "controlled_wing_phylo_results.rds")
    add_scalar("phylo_mean_wing_mm", r$mean_wing_mm, "controlled_wing_phylo_results.rds")
    add_scalar("phylo_decision_gate_scenario", r$decision_gate$scenario, "controlled_wing_phylo_results.rds")
    models_present <- if (!is.null(r$before_after)) paste(sort(unique(r$before_after$model)), collapse = ";") else ""
    add_scalar("phylo_models_present", models_present, "controlled_wing_phylo_results.rds")
  } else message("  (no controlled_wing_phylo_results.rds yet: figures fall back to the lme4 fast tier)")
  f <- out_path("controlled_wing_phylo_species_slopes.rds")
  if (file.exists(f)) {
    s <- readRDS(f)
    wcsv(bind_rows(s$M3, s$M3_cc, s$M0), "fig_species_slopes_phylo.csv")
    for (m in names(s$summary)) for (k in names(s$summary[[m]]))
      add_scalar(paste("slopes_phylo", m, k, sep = "."), s$summary[[m]][[k]], "controlled_wing_phylo_species_slopes.rds")
  } else message("  (no controlled_wing_phylo_species_slopes.rds yet: caterpillar figure falls back to the lme4 fast tier)")

  # --- brms Rubin-pooled results (Tiers 2-3, Bayesian cross-check), if the server run has been copied back ---
  f <- out_path("controlled_wing_brms_results.rds")
  if (file.exists(f)) {
    b <- readRDS(f)
    if (!is.null(b$species_slopes_M3)) wcsv(b$species_slopes_M3, "fig_species_slopes_brms_M3.csv")
    if (!is.null(b$before_after))      wcsv(b$before_after, "fig_controlled_before_after_brms.csv")
    add_scalar("brms_available", TRUE, "controlled_wing_brms_results.rds")
  } else message("  (no controlled_wing_brms_results.rds yet: brms Tier 2-3 cross-check not shown)")

  wcsv(bind_rows(scal), "fig_scalars.csv")
  invisible(fd)
}

fpath <- function(name, mode) file.path(FIG_DIR, paste0(name, "_", mode, ".png"))
save_fig <- function(name, mode, p, w = 7, h = 5) {
  tryCatch({
    ggsave(fpath(name, mode), p, width = w, height = h, dpi = 300, bg = "white")
    message("  saved ", fpath(name, mode))
  }, error = function(e) message("  FAILED ", name, " [", mode, "]: ", conditionMessage(e)))
}

# ============================================================================
# Figure builders (identical grammar to atlantic_drmsem.R Section 12, but fed
# from the saved data frames instead of live model objects)
# ============================================================================

# --- Annotated DAG (from standardised path table) --------------------------
# The wing-length response is drawn as TWO co-equal response nodes — mean(wing)
# and sd(wing) — because drmSEM models the mean and the residual SD as separate
# distributional responses (each with its own submodel). Mean-model paths point
# into mean(wing); the variance (sigma) paths point into sd(wing). Both node
# types share the "response" role colour and both receive ordinary directed
# arrows, so the residual-SD channel reads as a response on equal footing with
# the mean rather than as a secondary annotation on a single wing node.
fig_dag <- function(ps, mode) {
  pos <- data.frame(
    node  = c("scaled_yr","scaled_lat","Sex","scaled_tmean","arthro_obs_std","mean_wing","sd_wing"),
    x     = c(0, 0, 0, 1.2, 1.2, 2.6, 2.6),
    y     = c(3.2, 2.0, 0.6, 3.2, 1.3, 2.7, 0.9),
    role  = c("exogenous","exogenous","exogenous","mediator","mediator","response","response"),
    label = c("scaled_yr","scaled_lat","Sex","scaled_tmean","arthro_obs_std",
              "mean(wing)","sd(wing)"),
    stringsAsFactors = FALSE)
  xy <- function(n, col) pos[[col]][match(n, pos$node)]
  # Every path is a directed edge into a node. Retarget the two distributional
  # components of wing_length onto the two response nodes.
  e <- ps
  e$fromn <- ifelse(e$term == "SexMale", "Sex", e$from)
  e$ton   <- e$to
  e$ton[e$to == "wing_length" & e$component == "mu"]    <- "mean_wing"
  e$ton[e$to == "wing_length" & e$component == "sigma"] <- "sd_wing"
  e <- e[e$fromn %in% pos$node & e$ton %in% pos$node, ]
  e$x0 <- xy(e$fromn,"x"); e$y0 <- xy(e$fromn,"y")
  e$x1 <- xy(e$ton,"x");   e$y1 <- xy(e$ton,"y")
  e$linef <- ifelse(e$p.value < 0.05, "significant", "n.s.")
  # Stagger label position along the edge by target node so the crossing
  # Sex -> mean(wing) and year -> sd(wing) labels do not collide mid-plot.
  e$lf <- ifelse(e$ton == "sd_wing", 0.74, 0.56)
  present  <- unique(c(e$fromn, e$ton))
  pos_draw <- pos[pos$node %in% present, ]
  .clim <- max(0.6, stats::quantile(abs(e$std.estimate[e$fromn != "Sex"]), 0.95, na.rm = TRUE))
  ggplot() +
    geom_segment(data = e,
      aes(x = x0, y = y0, xend = x1, yend = y1, colour = std.estimate, linetype = linef),
      arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"), linewidth = 0.7) +
    geom_text(data = e,
      aes(x = x0 + lf*(x1-x0), y = y0 + lf*(y1-y0),
          label = sprintf("%.2f", std.estimate), colour = std.estimate),
      size = 2.7, fontface = "bold") +
    geom_label(data = pos_draw, aes(x = x, y = y, label = label, fill = role),
      colour = "black", size = 3.1) +
    scale_colour_gradient2(low = "#b2182b", mid = "grey75", high = "#1b7837",
      midpoint = 0, limits = c(-.clim, .clim), oob = scales::squish, name = "std. coef") +
    scale_fill_manual(values = c(exogenous = "#ECECEC", mediator = "#CDE7DD",
      response = "#FCE3C8"), name = NULL) +
    scale_linetype_manual(values = c(significant = "solid", `n.s.` = "dashed"), name = NULL) +
    coord_cartesian(xlim = c(-0.55, 3.2), ylim = c(0.1, 3.8), clip = "off") +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.box.margin = margin(t = 10, b = 15),
      legend.margin = margin(t = 5, b = 5),
      plot.margin = margin(t = 20, r = 40, b = 25, l = 40),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold", margin = margin(t = 5, b = 15))
    ) +
    labs(title = paste0("drmSEM path diagram [", mode,
                        "] — standardized coefficients; wing mean and SD are separate responses"))
}

# --- Coefficient forest (raw paths +/- 95% CI, faceted mu vs sigma) --------
fig_coef <- function(pd, mode) {
  pd$label <- paste(pd$from, "->", pd$to, "(", pd$term, ")")
  pd$lo  <- pd$estimate - 1.96 * pd$std.error
  pd$hi  <- pd$estimate + 1.96 * pd$std.error
  pd$sig <- ifelse(pd$p.value < 0.05, "p < 0.05", "n.s.")
  ggplot(pd, aes(estimate, reorder(label, estimate), colour = sig)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = lo, xmax = hi)) +
    facet_wrap(~ component, scales = "free", ncol = 1) +
    scale_colour_manual(values = c("p < 0.05" = "#1b7837", "n.s." = "grey60")) +
    labs(x = "coefficient (mu: mm per SD of predictor; sigma: log scale)",
         y = NULL, colour = NULL,
         title = paste0("drmSEM path coefficients [", mode, "]")) +
    theme_minimal()
}

# --- Effect-decomposition forest (year and temperature -> wing) ------------
fig_effects <- function(ey, et, mode) {
  mk <- function(e, lbl) { d <- as.data.frame(e); d$panel <- lbl; d }
  ed <- rbind(mk(ey, "scaled_yr -> wing"), mk(et, "scaled_tmean -> wing"))
  ed$quantity <- factor(ed$quantity,
    levels = c("total_path","direct","indirect","mean_mediated","distribution_mediated"))
  ggplot(ed, aes(estimate, quantity)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
    facet_wrap(~ panel, scales = "free_x") +
    labs(x = "effect on wing length (mm)", y = NULL,
         title = paste0("Effect decomposition [", mode, "]")) +
    theme_minimal()
}

# --- Distributional (sigma) channel: fold-change in residual SD ------------
fig_sigma <- function(pr, mode) {
  sr  <- pr[pr$component == "sigma", , drop = FALSE]
  lab <- c(scaled_yr = "Year (SD)", scaled_tmean = "Temperature (SD)")
  xx  <- seq(-2, 2, length.out = 60)
  cur <- do.call(rbind, lapply(seq_len(nrow(sr)), function(i) {
    b <- sr$estimate[i]; se <- sr$std.error[i]
    f_lo <- exp((b - 1.96 * se) * xx); f_hi <- exp((b + 1.96 * se) * xx)
    data.frame(
      predictor = ifelse(sr$term[i] %in% names(lab), lab[sr$term[i]], sr$term[i]),
      x = xx, fold = exp(b * xx), lo = pmin(f_lo, f_hi), hi = pmax(f_lo, f_hi))
  }))
  ggplot(cur, aes(x, fold)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
    geom_line(linewidth = 1, colour = "#762a83") +
    facet_wrap(~ predictor, scales = "free_x") +
    labs(x = "predictor (SD units)", y = "fold-change in residual SD of wing length",
         title = paste0("Distributional (sigma) channel [", mode, "]")) +
    theme_minimal()
}

# ============================================================================
# Drive: loop over requested modes
# ============================================================================
args  <- commandArgs(trailingOnly = TRUE)
if (!length(args) || "export" %in% args) {
  message("Exporting figure data for Manuscript/make_figures.py -> ", out_path("figure_data"))
  export_figure_data()
  if (identical(args, "export")) quit(save = "no", status = 0)
}
modes <- if (length(args)) setdiff(args, "export") else c("noarthro", "arthro")

for (mode in modes) {
  f <- file.path(OUTPUT_DIR, paste0("drmsem_results_", mode, ".rds"))
  if (!file.exists(f)) { message("skip [", mode, "]: ", basename(f), " not found"); next }
  res <- readRDS(f)
  ps  <- as.data.frame(res$paths_sdx)
  pr  <- as.data.frame(res$paths_raw)
  message("Figures for [", mode, "]  (N = ", res$n_records,
          ", generated ", res$generated, ")")
  save_fig("drmsem_dag",            mode, fig_dag(ps, mode),  9, 6.5)
  if (mode == "noarthro") {
    ms_dest <- file.path(.repo, "Manuscript", "images", "fig-drmsem.png")
    file.copy(fpath("drmsem_dag", mode), ms_dest, overwrite = TRUE)
    message("  copied drmsem_dag_noarthro.png -> Manuscript/images/fig-drmsem.png")
  }
  save_fig("drmsem_coef_forest",    mode, fig_coef(pr, mode), 7, 6)
  save_fig("drmsem_sigma_curves",   mode, fig_sigma(pr, mode), 8, 4)
  ey <- res$effects$yr; et <- res$effects$tmean
  if (!is.null(ey) && !is.null(et))
    save_fig("drmsem_effects_forest", mode, fig_effects(ey, et, mode), 9, 4)
  else message("  (no effect tables in ", basename(f), " — skipping effects forest)")
}
message("Done. Figures in ", FIG_DIR)
