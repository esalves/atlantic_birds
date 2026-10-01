# jirinec_benchmark.R
# ---------------------------------------------------------------------------
# WHY (internal review before the Oikos submission, 2026-09-30): the manuscript
# contrasts its trends with the central Amazonian study of Jirinec et al. (2021,
# Sci. Adv. 7: eabk1743) but never put the two on the same scale. This script
# derives the Amazonian benchmark from that paper's supplementary Table S6
# (data/external/jirinec2021_tableS6.csv) and plots it against our estimates.
#
# Benchmark. Table S6 gives, per species, model-estimated means in 1980 and 2019
# and their change in %. The community benchmark is the mean of the 77 species'
# changes divided by 3.9 decades (also for the 68 passerines). The paper reports
# warming at the site of +1.00 degC (wet season) and +1.65 degC (dry season)
# since 1966 (to 2019, 5.3 decades), i.e. 0.19-0.31 degC/decade; per-degC rates
# use that range. Ours: posterior means and 95% CrIs pooled over 20 trees
# (bayes_summary.rds): body mass (log scale, % per decade) and wing length
# (mm per decade, converted to % with the manuscript's effect-scale factor);
# per-degC rates divide by the regional warming rate (climate_trends.rds),
# without propagating its uncertainty.
#
# RUN:  Rscript Analysis/scripts/jirinec_benchmark.R
# OUT:  output/jirinec_benchmark.rds, Manuscript/images/fig-amazon.png
# ---------------------------------------------------------------------------
suppressMessages({ library(dplyr); library(ggplot2) })

.find_analysis_dir <- function() {
  d <- normalizePath(getwd(), winslash = "/")
  repeat {
    if (dir.exists(file.path(d, "data", "derived")) && dir.exists(file.path(d, "scripts"))) return(d)
    if (dir.exists(file.path(d, "Analysis", "data", "derived"))) return(normalizePath(file.path(d, "Analysis"), winslash = "/"))
    parent <- dirname(d); if (identical(parent, d)) break; d <- parent
  }
  stop("Could not locate Analysis/.")
}
ANALYSIS_DIR <- .find_analysis_dir()
out_path <- function(...) file.path(ANALYSIS_DIR, "output", ...)

# --- Amazonian benchmark ---------------------------------------------------------
jt <- read.csv(file.path(ANALYSIS_DIR, "data", "external", "jirinec2021_tableS6.csv"), comment.char = "#")
stopifnot(nrow(jt) == 77, sum(jt$mass_change_pct < 0) == 77, sum(jt$mass_sig == "decrease") == 36,
          sum(jt$wing_sig == "increase") == 22)
NONPASS <- c("Trochilidae", "Trogonidae", "Momotidae", "Galbulidae", "Bucconidae")
SPAN_DEC <- (2019 - 1980) / 10
WARM <- c(wet = 1.00, dry = 1.65) / ((2019 - 1966) / 10)          # degC per decade at the BDFFP
bench_of <- function(d) data.frame(
  trait = c("Body mass", "Wing length"),
  n_species = nrow(d),
  mean_pct_dec = c(mean(d$mass_change_pct), mean(d$wing_change_pct)) / SPAN_DEC,
  median_pct_dec = c(median(d$mass_change_pct), median(d$wing_change_pct)) / SPAN_DEC,
  min_pct_dec = c(min(d$mass_change_pct), min(d$wing_change_pct)) / SPAN_DEC,
  max_pct_dec = c(max(d$mass_change_pct), max(d$wing_change_pct)) / SPAN_DEC)
bench <- bind_rows(cbind(set = "all species", bench_of(jt)),
                   cbind(set = "passerines", bench_of(filter(jt, !family %in% NONPASS))))
bench$per_degC_lo <- bench$mean_pct_dec / WARM[["wet"]]          # slower warming -> larger per-degC rate
bench$per_degC_hi <- bench$mean_pct_dec / WARM[["dry"]]
print(bench, digits = 3)

# --- Atlantic Forest estimates (as quoted in the manuscript) ----------------------
bayes <- readRDS(out_path("bayes", "bayes_summary.rds"))
DEC   <- readRDS(out_path("bayes", "bayes_derived.rds"))$dec
es    <- readRDS(out_path("effect_scale.rds"))
cw    <- readRDS(out_path("controlled_wing_phylo_results.rds"))
clim  <- readRDS(out_path("climate_trends.rds"))
wpct  <- es$wing$using_scaling_sd["est", "pct_per_decade"] / cw$before_after$mm_per_decade[cw$before_after$model == "M0"][1]
warm_af <- clim$pooled$slope_per_decade[clim$pooled$variable == "tmean" & clim$pooled$weighting == "unweighted + year RE"]
row_of <- function(fit, scale) {
  r <- bayes$fixed[bayes$fixed$fit == fit & bayes$fixed$par == "b_scaled_yr", ]; stopifnot(nrow(r) == 1)
  f <- if (scale == "pct") function(x) 100 * (exp(x * DEC) - 1) else function(x) x * DEC * wpct
  c(est = f(r$estimate), lo = f(r$lower), hi = f(r$upper))
}
ours <- rbind(
  data.frame(trait = "Body mass", model = "unadjusted", t(row_of("atlantic_multitrait__001__y", "pct"))),
  data.frame(trait = "Body mass", model = "fully adjusted", t(row_of("atlantic_multitrait__003__y", "pct"))),
  data.frame(trait = "Wing length", model = "unadjusted", t(row_of("atlantic_parallel_controlled__001__conc.wing.length__M0", "mm"))),
  data.frame(trait = "Wing length", model = "fully adjusted", t(row_of("atlantic_parallel_controlled__011__conc.wing.length__M3", "mm"))))
print(ours, digits = 3)

saveRDS(list(generated = Sys.time(), source = "Jirinec et al. 2021 Sci. Adv. 7: eabk1743, Table S6",
             span_decades = SPAN_DEC, warming_bdffp_per_decade = WARM, warming_af_per_decade = warm_af,
             benchmark = bench, ours = ours, wing_pct_factor = wpct),
        out_path("jirinec_benchmark.rds"))

# --- figure ---------------------------------------------------------------------------
pb <- filter(bench, set == "passerines")
lev <- c("Central Amazon\n(68 passerines)", "Atlantic Forest,\nunadjusted", "Atlantic Forest,\nfully adjusted")
dA <- bind_rows(
  transmute(ours, trait, who = ifelse(model == "unadjusted", lev[2], lev[3]), est, lo, hi, panel = "A"),
  transmute(pb, trait, who = lev[1], est = mean_pct_dec, lo = min_pct_dec, hi = max_pct_dec, panel = "A"))
dB <- bind_rows(
  transmute(filter(ours, model == "fully adjusted"), trait, who = lev[3], est = est / warm_af, lo = lo / warm_af, hi = hi / warm_af, panel = "B"),
  transmute(pb, trait, who = lev[1], est = mean_pct_dec / mean(WARM), lo = pmin(per_degC_lo, per_degC_hi),
            hi = pmax(per_degC_lo, per_degC_hi), panel = "B"))
d <- bind_rows(dA, dB) %>% mutate(who = factor(who, levels = rev(lev)), trait = factor(trait, levels = c("Body mass", "Wing length")))
cols <- setNames(c("#C05621", "#A0AEC0", "#276749"), lev)
shp  <- setNames(c(18, 1, 16), lev)
xlab <- c(A = "Change (% per decade)", B = "Change (% per °C of local warming)")
mk <- function(p) {
  ggplot(filter(d, panel == p), aes(est, who, colour = who, shape = who)) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.8, orientation = "y") +
    geom_point(size = 3, fill = "white") +
    facet_wrap(~trait, ncol = 1) +
    scale_colour_manual(values = cols, guide = "none") + scale_shape_manual(values = shp, guide = "none") +
    scale_y_discrete(drop = FALSE) +                       # same rows in both panels
    labs(x = xlab[[p]], y = NULL, title = p) +
    theme_classic(base_size = 10, base_family = "Helvetica") +
    theme(plot.title = element_text(face = "bold", size = 11), strip.background = element_blank(),
          strip.text = element_text(face = "bold", hjust = 0), panel.spacing = unit(1, "lines"))
}
fig <- (function(a, b) { g <- gridExtra::arrangeGrob(a, b, ncol = 2, widths = c(1.15, 1)); g })(mk("A"), mk("B") + theme(axis.text.y = element_blank()))
fig_file <- file.path(dirname(ANALYSIS_DIR), "Manuscript", "images", "fig-amazon.png")
ggsave(fig_file, fig, width = 7.2, height = 4.0, dpi = 300, bg = "white")
message("wrote ", fig_file)
