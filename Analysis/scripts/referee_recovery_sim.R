# referee_recovery_sim.R
# ---------------------------------------------------------------------------
# WHY (internal review before the Oikos submission, 2026-09-30): capture year is
# strongly confounded with contributor and municipality, so the contributor and
# municipality intercepts of the fully adjusted model (M3) could absorb a real
# trend as easily as they remove an artefact. This script checks which model
# recovers a KNOWN trend when that trend is planted in the real sampling
# structure, and puts an interval on the "attenuation" of the wing decline.
#
# PART 1 - recovery simulation. The design (records, species, individuals,
# contributors, municipalities, years and covariates) is the real 8,478-record
# wing sample. A generating model G = M3 is fitted to the real data; new
# responses are built from G's fixed effects with the year slope replaced by a
# planted value (0, -0.5, -1.0 mm/decade), fresh species, individual and
# residual draws at G's variances, and contributor + municipality effects under
# three scenarios:
#   clean     fresh draws, independent of year: no artefact. Tests whether the
#             random intercepts absorb a real trend (bias of M3 towards zero).
#   real      the conditional modes (BLUPs) of G, so the real association of
#             contributor and municipality with year is kept.
#   artefact  fresh draws plus a contributor shift proportional to the
#             contributor's mean capture year, calibrated so that the artefact
#             alone gives M0 the observed M0 - M3 gap. No true trend unless planted.
# Each simulated data set is fitted with M0, M1, M3 and M5src (year within
# contributor), all lme4 REML with the same formulas as atlantic_parallel_controlled.R.
# Reported per scenario x planted value x model: mean estimate, bias, RMSE,
# coverage of the 95% Wald interval and the share of intervals excluding zero.
#
# PART 2 - attenuation interval. Cluster bootstrap over contributors on the real
# data (contributors resampled with replacement, all their records kept, copies
# relabelled): M0 and M3 refitted on each draw; attenuation = 1 - b(M3)/b(M0).
#
# The phylogenetic species intercept is omitted (lme4): it structures species
# intercepts, not the record-level year slope under test, and Tier 1 lme4 and the
# phylogenetic tiers agree on the year slope (atlantic_parallel_controlled.R).
#
# RUN (Totoro): Rscript Analysis/scripts/referee_recovery_sim.R --reps 200 --boot 200 --cores 20
#      smoke:   Rscript Analysis/scripts/referee_recovery_sim.R --smoke
# OUT: Analysis/output/referee_reruns/recovery_sim.rds (+ recovery_sim_summary.csv)
# ---------------------------------------------------------------------------

# --- data frames and model specs: sections 1-2 of atlantic_parallel_controlled.R ---
local({
  f <- c("Analysis/scripts/atlantic_parallel_controlled.R", "scripts/atlantic_parallel_controlled.R",
         "atlantic_parallel_controlled.R")
  f <- f[file.exists(f)][1]
  if (is.na(f)) stop("Run from the repo root, Analysis/ or Analysis/scripts/.")
  ex <- parse(f)
  stop_at <- which(vapply(ex, function(e) startsWith(deparse(e)[1], "YEAR_TERMS <-"), TRUE))[1]
  for (e in ex[seq_len(stop_at - 1)]) eval(e, envir = globalenv())
})
suppressMessages(library(parallel))

opt <- function(flag, default) { a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i)) default else if (i == length(a) || startsWith(a[i + 1], "--")) TRUE else a[i + 1] }
SMOKE <- isTRUE(opt("--smoke", FALSE))
REPS  <- if (SMOKE) 2L else as.integer(opt("--reps", 200))
BOOT  <- if (SMOKE) 2L else as.integer(opt("--boot", 200))
CORES <- if (SMOKE) 2L else as.integer(opt("--cores", max(1L, detectCores() - 1L)))
PLANTED <- c(0, -0.5, -1.0)                                   # mm per decade
OUT_DIR <- out_path("referee_reruns"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
tag <- if (SMOKE) "_SMOKE" else ""
RNGkind("L'Ecuyer-CMRG"); set.seed(20260930)

ctrl  <- lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
MMDEC <- 10 / SD_YEAR                                         # scaled_yr coefficient -> mm per decade
d     <- w_all
FITS  <- list(M0 = list(spec = SPECS$M0, term = "scaled_yr"),
              M1 = list(spec = SPECS$M1, term = "scaled_yr"),
              M3 = list(spec = SPECS$M3, term = "scaled_yr"),
              M5src = list(spec = SPECS$M5src, term = "yr_within_src"))

year_row <- function(fit, term) {
  co <- summary(fit)$coefficients
  if (!term %in% rownames(co)) return(c(est = NA_real_, se = NA_real_))
  c(est = co[term, 1] * MMDEC, se = co[term, 2] * MMDEC)
}
fit_all <- function(dat) {
  do.call(rbind, lapply(names(FITS), function(m) {
    fit <- tryCatch(suppressWarnings(lmer(lme4_formula(FITS[[m]]$spec), data = dat, REML = TRUE, control = ctrl)),
                    error = function(e) NULL)
    r <- if (is.null(fit)) c(est = NA_real_, se = NA_real_) else year_row(fit, FITS[[m]]$term)
    data.frame(model = m, est = r[["est"]], se = r[["se"]])
  }))
}

# --- generating model and its components ------------------------------------------
t0 <- Sys.time()
G  <- lmer(lme4_formula(SPECS$M3), data = d, REML = TRUE, control = ctrl)
M0 <- lmer(lme4_formula(SPECS$M0), data = d, REML = TRUE, control = ctrl)
vc <- as.data.frame(VarCorr(G))
sd_of <- function(grp, var = "(Intercept)") vc$sdcor[startsWith(vc$grp, grp) & vc$var1 == var & is.na(vc$var2)][1]
SD <- list(spp0 = sd_of("spp"), spp1 = sd_of("spp", "scaled_yr"), ind = sd_of("ind"),
           src = sd_of("src"), site = sd_of("site"), sigma = sigma(G))
X     <- getME(G, "X")
beta  <- fixef(G)
re_G  <- ranef(G)
blup_src  <- setNames(re_G$src[, 1],  rownames(re_G$src))
blup_site <- setNames(re_G$site[, 1], rownames(re_G$site))
obs <- list(M0 = unname(fixef(M0)["scaled_yr"]) * MMDEC, M3 = unname(beta["scaled_yr"]) * MMDEC)
gap_scaled <- unname(fixef(M0)["scaled_yr"] - beta["scaled_yr"])
# share of capture-year variance lying between contributors (record level)
r2_src <- var(d$yr_src_mean) / var(d$scaled_yr)
kappa  <- gap_scaled / r2_src      # artefact: src shift = kappa * contributor mean year
message(sprintf("[sim] observed M0 %.2f, M3 %.2f mm/decade; between-contributor share of year variance %.2f; kappa %.3f (%.0f s)",
                obs$M0, obs$M3, r2_src, kappa, as.numeric(Sys.time() - t0, units = "secs")))

idx <- list(spp = as.integer(d$spp), ind = as.integer(d$ind), src = as.integer(d$src), site = as.integer(d$site))
lev <- list(spp = nlevels(d$spp), ind = nlevels(d$ind), src = nlevels(d$src), site = nlevels(d$site))
yr_src_mean_lvl <- tapply(d$yr_src_mean, d$src, `[`, 1)

simulate_y <- function(planted_mm, scenario) {
  b <- beta; b["scaled_yr"] <- planted_mm / MMDEC
  mu <- drop(X %*% b)
  u0 <- rnorm(lev$spp, 0, SD$spp0); u1 <- rnorm(lev$spp, 0, SD$spp1)
  mu <- mu + u0[idx$spp] + u1[idx$spp] * d$scaled_yr + rnorm(lev$ind, 0, SD$ind)[idx$ind]
  if (scenario == "real") {
    u_src <- blup_src[levels(d$src)]; u_site <- blup_site[levels(d$site)]
  } else {
    u_src <- rnorm(lev$src, 0, SD$src); u_site <- rnorm(lev$site, 0, SD$site)
    if (scenario == "artefact") u_src <- u_src + kappa * yr_src_mean_lvl
  }
  mu + u_src[idx$src] + u_site[idx$site] + rnorm(nrow(d), 0, SD$sigma)
}

cells <- expand.grid(scenario = c("clean", "real", "artefact"), planted = PLANTED, rep = seq_len(REPS),
                     stringsAsFactors = FALSE)
message(sprintf("[sim] %d simulated data sets x %d models on %d cores", nrow(cells), length(FITS), CORES))
sim_rows <- mclapply(seq_len(nrow(cells)), function(k) {
  cl <- cells[k, ]; dk <- d; dk$conc.wing.length <- simulate_y(cl$planted, cl$scenario)
  cbind(cl, fit_all(dk))
}, mc.cores = CORES, mc.set.seed = TRUE)
bad <- vapply(sim_rows, inherits, TRUE, "try-error")
if (any(bad)) warning(sum(bad), " simulation cells failed and are dropped")
sim <- do.call(rbind, sim_rows[!bad])
sim$lo <- sim$est - 1.96 * sim$se; sim$hi <- sim$est + 1.96 * sim$se

summary_tbl <- sim %>% filter(is.finite(est)) %>%
  group_by(scenario, planted, model) %>%
  summarise(n = n(), mean_est = mean(est), bias = mean(est - planted), rmse = sqrt(mean((est - planted)^2)),
            coverage = mean(lo <= planted & hi >= planted), excl_zero = mean(hi < 0 | lo > 0),
            mean_se = mean(se), .groups = "drop") %>%
  arrange(scenario, planted, match(model, names(FITS)))
print(as.data.frame(summary_tbl), digits = 3)

# --- PART 2: cluster bootstrap of the attenuation ----------------------------------
src_lv <- levels(d$src)
boot_rows <- mclapply(seq_len(BOOT), function(b) {
  pick <- sample(src_lv, length(src_lv), replace = TRUE)
  db <- bind_rows(lapply(seq_along(pick), function(k) {
    x <- d[d$src == pick[k], ]; x$src <- paste0(pick[k], "#", k); x$ind <- paste0(x$ind, "#", k); x }))
  db <- db %>% mutate(src = factor(src), ind = factor(ind)) %>% droplevels()
  f0 <- tryCatch(suppressWarnings(lmer(lme4_formula(SPECS$M0), data = db, REML = TRUE, control = ctrl)), error = function(e) NULL)
  f3 <- tryCatch(suppressWarnings(lmer(lme4_formula(SPECS$M3), data = db, REML = TRUE, control = ctrl)), error = function(e) NULL)
  if (is.null(f0) || is.null(f3)) return(data.frame(b = b, M0 = NA_real_, M3 = NA_real_))
  data.frame(b = b, M0 = fixef(f0)[["scaled_yr"]] * MMDEC, M3 = fixef(f3)[["scaled_yr"]] * MMDEC)
}, mc.cores = CORES, mc.set.seed = TRUE)
boot <- do.call(rbind, boot_rows[!vapply(boot_rows, inherits, TRUE, "try-error")])
boot$attenuation <- 1 - boot$M3 / boot$M0
ok <- is.finite(boot$attenuation)
atten <- list(observed = 1 - obs$M3 / obs$M0, n_boot = sum(ok),
              median = median(boot$attenuation[ok]),
              ci95 = unname(quantile(boot$attenuation[ok], c(.025, .975))),
              diff_ci95 = unname(quantile((boot$M3 - boot$M0)[ok], c(.025, .975))))
message(sprintf("[boot] attenuation observed %.0f%%, bootstrap median %.0f%%, 95%% interval %.0f%% to %.0f%% (%d draws)",
                100 * atten$observed, 100 * atten$median, 100 * atten$ci95[1], 100 * atten$ci95[2], atten$n_boot))

res <- list(generated = Sys.time(), smoke = SMOKE, reps = REPS, boot = BOOT, planted_mm_per_decade = PLANTED,
            generating_model = "M3 (lme4 REML) on w_all", observed = obs, sd = SD,
            artefact = list(kappa_scaled = kappa, r2_between_contributor = r2_src, gap_mm = gap_scaled * MMDEC),
            summary = as.data.frame(summary_tbl), sims = sim, attenuation = atten, bootstrap = boot,
            secs = as.numeric(Sys.time() - t0, units = "secs"), lme4_version = as.character(packageVersion("lme4")))
saveRDS(res, file.path(OUT_DIR, paste0("recovery_sim", tag, ".rds")))
write.csv(summary_tbl, file.path(OUT_DIR, paste0("recovery_sim_summary", tag, ".csv")), row.names = FALSE)
message(sprintf("[done] %.1f min -> %s", res$secs / 60, file.path(OUT_DIR, paste0("recovery_sim", tag, ".rds"))))
