# make_figures.R
# ---------------------------------------------------------------------------
# Regenerate all drmSEM result figures from the saved Analysis/output/*.rds
# WITHOUT refitting any model. Reads drmsem_results_<mode>.rds (which store
# paths_raw, paths_sdx and the effect tables as plain data frames) and redraws
# the four figures per mode into Analysis/figures/.
#
# Only needs ggplot2 (+ scales, grid — both ggplot2 deps / base). No drmTMB,
# drmSEM, or fitted model objects required, so it runs fast anywhere R + ggplot2
# are installed.
#
#   Rscript make_figures.R            # both modes (noarthro, arthro)
#   Rscript make_figures.R arthro     # one mode only
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(ggplot2))

# --- locate repo (Analysis/output + Analysis/figures) ----------------------
.repo <- local({
  d <- getwd()
  for (i in 0:10) {
    if (file.exists(file.path(d, "atlantic_birds.Rproj")) ||
        dir.exists(file.path(d, "Analysis", "output"))) return(d)
    d <- dirname(d)
  }
  stop("Could not find the repo root (looked for atlantic_birds.Rproj / Analysis/output) from ", getwd())
})
OUTPUT_DIR <- file.path(.repo, "Analysis", "output")
FIG_DIR    <- file.path(.repo, "Analysis", "figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

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
modes <- if (length(args)) args else c("noarthro", "arthro")

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
