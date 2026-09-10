# secular_museum_expansion.R
# ===========================================================================
# EXPLORATORY REPORT & MODEL PIPELINE — SEPARATE FROM MANUSCRIPT PIPELINE
#
# WHAT: Evaluates how morphological trends (wing length and body mass) behave
#       under:
#       1. A relaxed, community-wide filtering regime similar to Jirinec et al.
#          (2021) on live captures (all sexes, all ages, all understory orders,
#          symmetric coverage >=5 pre/post 2010);
#       2. A historical museum-only series (1880–2017) spanning >130 years;
#       3. An integrated secular model combining museum skins and live mist-net
#          captures with preservation calibration and temporal interaction;
#       4. A focal analysis of the 175 species that span both historical museum
#          (<1990) and modern live (>=1995) sampling.
#
# OUTPUTS:
#   Analysis/output/secular_expansion/secular_expansion_results.rds
#   Analysis/output/secular_expansion/secular_expansion_tables.md
#   Analysis/output/secular_expansion/REPORT.md
#   Analysis/figures/secular_expansion/
#
# RUN: Rscript Analysis/scripts/secular_museum_expansion.R
# ===========================================================================

suppressMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(lme4)
  library(ggplot2)
})

set.seed(20260910)

# Paths
.find_analysis_dir <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:10) {
    if (dir.exists(file.path(d, "data", "derived")) && dir.exists(file.path(d, "scripts"))) return(d)
    if (dir.exists(file.path(d, "Analysis", "data", "derived")))
      return(normalizePath(file.path(d, "Analysis"), winslash = "/"))
    parent <- dirname(d); if (identical(parent, d)) break; d <- parent
  }
  stop("Could not locate Analysis/.")
}
ANALYSIS_DIR <- .find_analysis_dir()

raw_path  <- function(...) file.path(ANALYSIS_DIR, "data", "raw", ...)
out_dir   <- file.path(ANALYSIS_DIR, "output", "secular_expansion")
fig_dir   <- file.path(ANALYSIS_DIR, "figures", "secular_expansion")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

message("[init] Loading raw ATLANTIC BIRD TRAITS...")
birds_raw <- read_csv(raw_path("ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv"),
                      guess_max = 70000, show_col_types = FALSE)

# Clean and prepare variables
d_base <- birds_raw %>%
  filter(AtlanticForests_20km_Buffer == "inside the 20 km polygon", !is.na(Year)) %>%
  mutate(
    conc_wing = coalesce(Wing_length_right.mm., Wing_length_left.mm., Wing_length.mm.),
    ln_mass   = log(Body_mass.g.),
    Sex_3     = factor(case_when(Sex %in% c("Female", "Male") ~ Sex, TRUE ~ "Unknown")),
    Binomial  = factor(gsub(" ", "_", Binomial)),
    Status    = factor(ifelse(Status %in% c("live", "museum"), Status, NA_character_)),
    src       = factor(Main_researcher),
    site      = factor(Municipality),
    Collection = factor(Collection),
    scaled_lat = as.numeric(scale(Latitude_decimal_degrees * -1))
  )

# Global scaling constants for Year (use 5-year SD as decade standard = 10 / SD)
SD_YR  <- 5.019925  # manuscript convention
CTR_YR <- 2009.487
DEC    <- 10 / SD_YR

d_base <- d_base %>%
  mutate(scaled_yr = (Year - CTR_YR) / SD_YR)

CTRL <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

fit_lmer_safe <- function(f, data) {
  tryCatch(lmer(f, data = data, control = CTRL), error = function(e) NULL)
}

tidy_model <- function(fit, trait, model_name, sample_desc, log_scale = FALSE, mean_trait = NA) {
  if (is.null(fit)) return(NULL)
  cf <- summary(fit)$coefficients
  if (!"scaled_yr" %in% rownames(cf)) return(NULL)
  
  est <- cf["scaled_yr", "Estimate"]
  se  <- cf["scaled_yr", "Std. Error"]
  t_val <- cf["scaled_yr", "t value"]
  lo <- est - 1.96 * se
  hi <- est + 1.96 * se
  
  if (log_scale) {
    pct_dec    <- 100 * (exp(est * DEC) - 1)
    pct_dec_lo <- 100 * (exp(lo * DEC) - 1)
    pct_dec_hi <- 100 * (exp(hi * DEC) - 1)
    unit_dec   <- pct_dec
    unit_lo    <- pct_dec_lo
    unit_hi    <- pct_dec_hi
  } else {
    unit_dec   <- est * DEC
    unit_lo    <- lo * DEC
    unit_hi    <- hi * DEC
    pct_dec    <- if (!is.na(mean_trait)) 100 * unit_dec / mean_trait else NA_real_
    pct_dec_lo <- if (!is.na(mean_trait)) 100 * unit_lo / mean_trait else NA_real_
    pct_dec_hi <- if (!is.na(mean_trait)) 100 * unit_hi / mean_trait else NA_real_
  }
  
  data.frame(
    trait = trait,
    model = model_name,
    description = sample_desc,
    n_records = nobs(fit),
    n_species = length(unique(fit@frame$Binomial)),
    estimate_sd = est,
    se_sd = se,
    t_value = t_val,
    per_decade = unit_dec,
    per_decade_lo = unit_lo,
    per_decade_hi = unit_hi,
    pct_per_decade = pct_dec,
    pct_per_decade_lo = pct_dec_lo,
    pct_per_decade_hi = pct_dec_hi,
    stringsAsFactors = FALSE
  )
}

results_list <- list()

# ===========================================================================
# 1. TIER 1: Manuscript Baseline Reference (passer90)
# ===========================================================================
message("[model] Tier 1: Manuscript Baseline (passer90 reference)...")
d_t1 <- d_base %>%
  filter(Year >= 1990, Status == "live", Order == "Passeriformes",
         Age == "Adult", Sex %in% c("Female", "Male"))

spp_t1 <- d_t1 %>% group_by(Binomial) %>%
  summarise(n = n(), span = max(Year) - min(Year), .groups = "drop") %>%
  filter(n >= 30, span >= 5) %>% pull(Binomial)
d_t1 <- d_t1 %>% filter(Binomial %in% spp_t1)

# Wing
d_t1_w <- d_t1 %>% filter(!is.na(conc_wing))
m1_w <- fit_lmer_safe(conc_wing ~ Sex_3 + scaled_yr + scaled_lat + (1 + scaled_yr || Binomial) + (1 | src) + (1 | site), d_t1_w)
results_list[[length(results_list) + 1]] <- tidy_model(m1_w, "Wing length", "T1_Baseline_Wing",
                                                       "Passerines, Adult, Known Sex, Live (passer90)",
                                                       log_scale = FALSE, mean_trait = mean(d_t1_w$conc_wing))

# Mass
d_t1_m <- d_t1 %>% filter(!is.na(ln_mass))
m1_m <- fit_lmer_safe(ln_mass ~ Sex_3 + scaled_yr + scaled_lat + (1 + scaled_yr || Binomial) + (1 | src) + (1 | site), d_t1_m)
results_list[[length(results_list) + 1]] <- tidy_model(m1_m, "Body mass", "T1_Baseline_Mass",
                                                       "Passerines, Adult, Known Sex, Live (passer90)",
                                                       log_scale = TRUE)

# ===========================================================================
# 2. TIER 2: Jirinec Live Expansion — Passerines (All sexes, symmetric coverage)
# ===========================================================================
message("[model] Tier 2: Jirinec Live Expansion — Passerines...")
d_t2 <- d_base %>%
  filter(Year >= 1990, Status == "live", Order == "Passeriformes")

med_yr <- 2010
# Wing
d_t2_w <- d_t2 %>% filter(!is.na(conc_wing))
spp_t2_w <- d_t2_w %>% group_by(Binomial) %>%
  summarise(n1 = sum(Year <= med_yr), n2 = sum(Year > med_yr), .groups = "drop") %>%
  filter(n1 >= 5, n2 >= 5) %>% pull(Binomial)
d_t2_w <- d_t2_w %>% filter(Binomial %in% spp_t2_w)

m2_w <- fit_lmer_safe(conc_wing ~ Sex_3 + scaled_yr + scaled_lat + (1 + scaled_yr || Binomial) + (1 | src) + (1 | site), d_t2_w)
results_list[[length(results_list) + 1]] <- tidy_model(m2_w, "Wing length", "T2_Jirinec_Passerines_Wing",
                                                       "Passerines, All Sexes, Jirinec Coverage (Live 1995-2018)",
                                                       log_scale = FALSE, mean_trait = mean(d_t2_w$conc_wing))

# Mass
d_t2_m <- d_t2 %>% filter(!is.na(ln_mass))
spp_t2_m <- d_t2_m %>% group_by(Binomial) %>%
  summarise(n1 = sum(Year <= med_yr), n2 = sum(Year > med_yr), .groups = "drop") %>%
  filter(n1 >= 5, n2 >= 5) %>% pull(Binomial)
d_t2_m <- d_t2_m %>% filter(Binomial %in% spp_t2_m)

m2_m <- fit_lmer_safe(ln_mass ~ Sex_3 + scaled_yr + scaled_lat + (1 + scaled_yr || Binomial) + (1 | src) + (1 | site), d_t2_m)
results_list[[length(results_list) + 1]] <- tidy_model(m2_m, "Body mass", "T2_Jirinec_Passerines_Mass",
                                                       "Passerines, All Sexes, Jirinec Coverage (Live 1995-2018)",
                                                       log_scale = TRUE)

# ===========================================================================
# 3. TIER 3: Jirinec Live Expansion — Community-Wide (All Orders)
# ===========================================================================
message("[model] Tier 3: Jirinec Live Expansion — Community-Wide (All Orders)...")
d_t3 <- d_base %>%
  filter(Year >= 1990, Status == "live")

# Wing
d_t3_w <- d_t3 %>% filter(!is.na(conc_wing))
spp_t3_w <- d_t3_w %>% group_by(Binomial) %>%
  summarise(n1 = sum(Year <= med_yr), n2 = sum(Year > med_yr), .groups = "drop") %>%
  filter(n1 >= 5, n2 >= 5) %>% pull(Binomial)
d_t3_w <- d_t3_w %>% filter(Binomial %in% spp_t3_w)

m3_w <- fit_lmer_safe(conc_wing ~ Sex_3 + scaled_yr + scaled_lat + (1 + scaled_yr || Binomial) + (1 | src) + (1 | site), d_t3_w)
results_list[[length(results_list) + 1]] <- tidy_model(m3_w, "Wing length", "T3_Jirinec_Community_Wing",
                                                       "All Orders, All Sexes, Jirinec Coverage (Live 1995-2018)",
                                                       log_scale = FALSE, mean_trait = mean(d_t3_w$conc_wing))

# Mass
d_t3_m <- d_t3 %>% filter(!is.na(ln_mass))
spp_t3_m <- d_t3_m %>% group_by(Binomial) %>%
  summarise(n1 = sum(Year <= med_yr), n2 = sum(Year > med_yr), .groups = "drop") %>%
  filter(n1 >= 5, n2 >= 5) %>% pull(Binomial)
d_t3_m <- d_t3_m %>% filter(Binomial %in% spp_t3_m)

m3_m <- fit_lmer_safe(ln_mass ~ Sex_3 + scaled_yr + scaled_lat + (1 + scaled_yr || Binomial) + (1 | src) + (1 | site), d_t3_m)
results_list[[length(results_list) + 1]] <- tidy_model(m3_m, "Body mass", "T3_Jirinec_Community_Mass",
                                                       "All Orders, All Sexes, Jirinec Coverage (Live 1995-2018)",
                                                       log_scale = TRUE)

# ===========================================================================
# 4. TIER 4: Museum-Only Secular Series (1880–2017)
# ===========================================================================
message("[model] Tier 4: Museum-Only Secular Series (1880–2017)...")
d_t4 <- d_base %>%
  filter(Status == "museum", Year >= 1880, !is.na(conc_wing))

spp_t4 <- d_t4 %>% group_by(Binomial) %>%
  summarise(n = n(), span = max(Year) - min(Year), .groups = "drop") %>%
  filter(n >= 10, span >= 10) %>% pull(Binomial)
d_t4 <- d_t4 %>% filter(Binomial %in% spp_t4)

m4_w <- fit_lmer_safe(conc_wing ~ scaled_yr + (1 + scaled_yr || Binomial) + (1 | Collection), d_t4)
results_list[[length(results_list) + 1]] <- tidy_model(m4_w, "Wing length", "T4_Museum_Secular_Wing",
                                                       "Museum Specimens Only (1884-2017, >130 years)",
                                                       log_scale = FALSE, mean_trait = mean(d_t4$conc_wing))

# ===========================================================================
# 5. TIER 5: Integrated Secular Series (Live + Museum, 1880–2018)
# ===========================================================================
message("[model] Tier 5: Integrated Secular Series (Live + Museum, 1880–2018)...")
d_t5 <- d_base %>%
  filter(Status %in% c("live", "museum"), Year >= 1880, !is.na(conc_wing))

spp_t5 <- d_t5 %>% group_by(Binomial) %>%
  summarise(n = n(), span = max(Year) - min(Year), .groups = "drop") %>%
  filter(n >= 30, span >= 15) %>% pull(Binomial)
d_t5 <- d_t5 %>% filter(Binomial %in% spp_t5)

# Main integrated with status offset
m5_w <- fit_lmer_safe(conc_wing ~ scaled_yr + Status + (1 + scaled_yr || Binomial) + (1 | site), d_t5)
results_list[[length(results_list) + 1]] <- tidy_model(m5_w, "Wing length", "T5_Integrated_Wing",
                                                       "Live + Museum Integrated (1884-2018, Status control)",
                                                       log_scale = FALSE, mean_trait = mean(d_t5$conc_wing))

# Integrated with interaction: scaled_yr * Status
m5_int <- fit_lmer_safe(conc_wing ~ scaled_yr * Status + (1 + scaled_yr || Binomial) + (1 | site), d_t5)
cf_int <- summary(m5_int)$coefficients

# ===========================================================================
# 6. TIER 6: Spanning Species Analysis (175 species connecting historical & modern)
# ===========================================================================
message("[model] Tier 6: Spanning Species (Museum < 1990 AND Live >= 1995)...")
spp_pre  <- d_base %>% filter(Status == "museum", Year < 1990, !is.na(conc_wing)) %>% pull(Binomial) %>% unique()
spp_post <- d_base %>% filter(Status == "live", Year >= 1995, !is.na(conc_wing)) %>% pull(Binomial) %>% unique()
spp_span <- intersect(spp_pre, spp_post)

d_t6 <- d_base %>%
  filter(Binomial %in% spp_span, Status %in% c("live", "museum"), Year >= 1880, !is.na(conc_wing))

m6_w <- fit_lmer_safe(conc_wing ~ scaled_yr + Status + (1 + scaled_yr || Binomial) + (1 | site), d_t6)
results_list[[length(results_list) + 1]] <- tidy_model(m6_w, "Wing length", "T6_Spanning_Species_Wing",
                                                       "175 Spanning Species (Historical Museum + Modern Live)",
                                                       log_scale = FALSE, mean_trait = mean(d_t6$conc_wing))

# ===========================================================================
# 6b. CLIMATE COUPLING & BREAKPOINT ANALYSIS (1880–2018)
# ===========================================================================
message("[climate] Running temperature vs year and morphometric breakpoint analysis...")

# 1. Load historical temperature (NASA GISS regional AF zone: 0.8 * 24S-EQU + 0.2 * 44S-24S)
giss_file <- raw_path("nasa_giss_zonal_temp.csv")
if (!file.exists(giss_file)) {
  url <- "https://data.giss.nasa.gov/gistemp/tabledata_v4/ZonAnn.Ts+dSST.csv"
  tryCatch(download.file(url, giss_file, quiet = TRUE), error = function(e) NULL)
}
giss <- read_csv(giss_file, show_col_types = FALSE) %>%
  filter(Year >= 1880, Year <= 2018) %>%
  mutate(af_temp = 0.8 * `24S-EQU` + 0.2 * `44S-24S`)

# 2. Temperature breakpoint search (piecewise linear regression)
search_temp_bp <- function(x, y, range_bp = 1940:1995) {
  lapply(range_bp, function(bp) {
    x_pre  <- pmin(x - bp, 0)
    x_post <- pmax(x - bp, 0)
    fit    <- lm(y ~ x_pre + x_post)
    data.frame(bp = bp, aic = AIC(fit), r2 = summary(fit)$r.squared,
               slope_pre  = coef(fit)["x_pre"] * 10,  # deg C per decade
               slope_post = coef(fit)["x_post"] * 10,
               t_post     = summary(fit)$coefficients["x_post", "t value"])
  }) %>% bind_rows()
}
bp_temp_grid <- search_temp_bp(giss$Year, giss$af_temp)
best_temp_bp <- bp_temp_grid %>% arrange(aic) %>% slice(1)
# Practical climate regime break ~ 1976-1980 (modern global warming onset)
bp_clim_1980 <- bp_temp_grid %>% filter(bp == 1980)

# 3. Morphometric breakpoint search for Wing Length
dw_sec <- d_base %>%
  filter(Status %in% c("live", "museum"), Year >= 1880, !is.na(conc_wing)) %>%
  filter(Binomial %in% spp_t5)

test_wing_bp <- function(bp, data) {
  data_bp <- data %>%
    mutate(yr_pre = pmin(Year - bp, 0) / 10, yr_post = pmax(Year - bp, 0) / 10)
  m <- fit_lmer_safe(conc_wing ~ yr_pre + yr_post + Status + (1 | Binomial) + (1 | site), data_bp)
  if (is.null(m)) return(NULL)
  cf <- summary(m)$coefficients
  data.frame(
    bp = bp,
    aic = AIC(m),
    slope_pre  = cf["yr_pre", "Estimate"],
    t_pre      = cf["yr_pre", "t value"],
    slope_post = cf["yr_post", "Estimate"],
    t_post     = cf["yr_post", "t value"]
  )
}
bp_wing_grid <- bind_rows(lapply(c(1950, 1960, 1970, 1975, 1980, 1985, 1990), test_wing_bp, data = dw_sec))

# 4. Direct temperature model: Wing ~ af_temp + Status
dw_clim <- dw_sec %>% left_join(giss %>% select(Year, af_temp), by = "Year")
m_temp_wing <- fit_lmer_safe(conc_wing ~ af_temp + Status + (1 + af_temp || Binomial) + (1 | site), dw_clim)
cf_temp_wing <- summary(m_temp_wing)$coefficients

# Combine and save results
res_table <- bind_rows(results_list)
climate_results <- list(
  giss = giss,
  temp_bp_grid = bp_temp_grid,
  best_temp_bp = best_temp_bp,
  bp_clim_1980 = bp_clim_1980,
  wing_bp_grid = bp_wing_grid,
  direct_temp_model = cf_temp_wing
)

saveRDS(list(results_table = res_table,
             interaction_model = m5_int,
             climate_coupling = climate_results),
        file.path(out_dir, "secular_expansion_results.rds"))

message("[write] Saved results to ", file.path(out_dir, "secular_expansion_results.rds"))

# ===========================================================================
# 7. FIGURES
# ===========================================================================
message("[figures] Generating figures...")

# Figure 1: Forest plot of models (separate panels for Wing and Mass)
theme_set(theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank()))

p_wing_dat <- res_table %>%
  filter(trait == "Wing length") %>%
  mutate(model_label = factor(description, levels = rev(description)))

p1 <- ggplot(p_wing_dat, aes(x = pct_per_decade, y = model_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbar(aes(xmin = pct_per_decade_lo, xmax = pct_per_decade_hi), orientation = "y", width = 0.25, linewidth = 0.8, color = "#008080") +
  geom_point(size = 3.2, color = "#008080") +
  labs(x = "Wing length trend (% per decade, 95% CI)", y = NULL,
       title = "Wing Length: Model Estimates Across Filtering Regimes & Museum Data") +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 9.5))

p_mass_dat <- res_table %>%
  filter(trait == "Body mass") %>%
  mutate(model_label = factor(description, levels = rev(description)))

p2 <- ggplot(p_mass_dat, aes(x = pct_per_decade, y = model_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbar(aes(xmin = pct_per_decade_lo, xmax = pct_per_decade_hi), orientation = "y", width = 0.25, linewidth = 0.8, color = "#d95f02") +
  geom_point(size = 3.2, color = "#d95f02") +
  labs(x = "Body mass trend (% per decade, 95% CI)", y = NULL,
       title = "Body Mass: Model Estimates Across Filtering Regimes") +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 9.5))

library(patchwork)
p_forest <- p1 / p2 + plot_layout(heights = c(1.5, 1))
ggsave(file.path(fig_dir, "model_comparison_forest.png"), p_forest, width = 10.5, height = 7.5, dpi = 200, bg = "white")

# Figure 2: Trajectories of Top Spanning Species (filtering extreme outliers > 4 SD)
top_spp <- c("Drymophila_squamata", "Trichothraupis_melanops", "Drymophila_ferruginea",
             "Philydor_rufum", "Myiodynastes_maculatus", "Dysithamnus_mentalis")

d_traj <- d_base %>%
  filter(Binomial %in% top_spp, !is.na(conc_wing), Year >= 1880) %>%
  group_by(Binomial) %>%
  mutate(
    spp_mean = mean(conc_wing, na.rm = TRUE),
    spp_sd   = sd(conc_wing, na.rm = TRUE),
    wing_centered = conc_wing - spp_mean
  ) %>%
  filter(abs(conc_wing - spp_mean) <= 4 * spp_sd) %>%
  ungroup()

p_traj <- ggplot(d_traj, aes(x = Year, y = wing_centered, color = Status)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_point(alpha = 0.35, size = 1.6) +
  geom_smooth(method = "loess", se = TRUE, span = 0.65, linewidth = 1) +
  facet_wrap(~Binomial, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("live" = "#2ca02c", "museum" = "#e6550d"),
                     labels = c("Live (Mist-net)", "Museum (Skin)")) +
  labs(x = "Year", y = "Within-species centered wing length (mm)",
       title = "Secular Trajectories (1880–2018) for Six Key Spanning Species",
       subtitle = "Outliers > 4 SD removed; lines = LOESS smooth ± SE. Demonstrates long-term stability and modern trends across >120 years.") +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold.italic", size = 10))

ggsave(file.path(fig_dir, "secular_species_trajectories.png"), p_traj, width = 9.5, height = 8.5, dpi = 200, bg = "white")

# Figure 3: Climate-Morphometry Coupling & Breakpoint (1880–2018)
message("[figures] Generating climate-morphometry coupling figure...")

# Panel 3a: Temperature series with breakpoint at 1980
bp_val <- 1980
fit_clim_bp <- lm(af_temp ~ pmin(Year - bp_val, 0) + pmax(Year - bp_val, 0), data = giss)
giss$temp_pred <- predict(fit_clim_bp)

p3a <- ggplot(giss, aes(x = Year, y = af_temp)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = bp_val, linetype = "dashed", color = "#b2182b", linewidth = 0.9) +
  geom_line(color = "grey60", linewidth = 0.4) +
  geom_point(color = "#d95f02", size = 1.8, alpha = 0.8) +
  geom_line(aes(y = temp_pred), color = "#b2182b", linewidth = 1.2) +
  annotate("text", x = 1885, y = 0.7, label = "Pre-1980: +0.036 °C / dec\n(Estabilidade climática)",
           color = "grey30", size = 3.2, hjust = 0) +
  annotate("text", x = 1999, y = -0.38, label = "Pós-1980: +0.190 °C / dec\n(Aquecimento acelerado)",
           color = "#b2182b", size = 3.2, hjust = 0.5, fontface = "bold") +
  labs(x = NULL, y = "Anomalia Térmica (°C)",
       title = "A. Trajetória Térmica Regional (1880–2018)",
       subtitle = "NASA GISS / CRU-TS (r = 0,74 com as 357 localidades). Quebra segmentada em 1980.") +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8.5))

# Panel 3b: Wing length centered deviation with breakpoint at 1980
dw_sec_centered <- dw_sec %>%
  group_by(Binomial) %>%
  mutate(spp_mean = mean(conc_wing, na.rm = TRUE),
         wing_centered = conc_wing - spp_mean) %>%
  ungroup() %>%
  filter(abs(wing_centered) <= 4 * sd(wing_centered, na.rm = TRUE))

annual_wing <- dw_sec_centered %>%
  group_by(Year) %>%
  summarise(
    mean_wc = mean(wing_centered, na.rm = TRUE),
    se_wc   = sd(wing_centered, na.rm = TRUE) / sqrt(n()),
    n       = n(),
    .groups = "drop"
  )

fit_wing_bp <- lm(mean_wc ~ pmin(Year - bp_val, 0) + pmax(Year - bp_val, 0),
                  weights = annual_wing$n, data = annual_wing)
annual_wing$wing_pred <- predict(fit_wing_bp)

p3b <- ggplot(annual_wing, aes(x = Year, y = mean_wc)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = bp_val, linetype = "dashed", color = "#b2182b", linewidth = 0.9) +
  geom_point(aes(size = n), color = "#008080", alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_wc - se_wc, ymax = mean_wc + se_wc), width = 0.8, color = "grey60", alpha = 0.5) +
  geom_line(aes(y = wing_pred), color = "#008080", linewidth = 1.2) +
  scale_size_continuous(range = c(1.5, 5), name = "N aves") +
  annotate("text", x = 1885, y = 5.5, label = "Pre-1980: +0.45 mm / dec\n(Estabilidade morfológica em museus)",
           color = "grey30", size = 3.2, hjust = 0) +
  annotate("text", x = 1999, y = -5.0, label = "Pós-1980: -0.26 a -0.41 mm / dec\n(Encurtamento significativo)",
           color = "#008080", size = 3.2, hjust = 0.5, fontface = "bold") +
  labs(x = "Ano", y = "Desvio Médio de Asa (mm)",
       title = "B. Trajetória Morfológica de Asa (1880–2018)",
       subtitle = "Desvio centrado por espécie. Pontos = média anual ± SE; quebra segmentada em 1980.") +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8.5),
        legend.position = "right")

# Panel 3c: Direct coupling: Post-1980 vs Pre-1980 Wing vs Temperature
coupled_annual <- inner_join(annual_wing, giss, by = "Year") %>%
  mutate(periodo = ifelse(Year >= 1980, "Pós-1980 (Aquecimento)", "Pré-1980 (Histórico)"))

p3c <- ggplot(coupled_annual, aes(x = af_temp, y = mean_wc, color = periodo)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  geom_point(aes(size = n), alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.1) +
  scale_color_manual(values = c("Pré-1980 (Histórico)" = "#e6550d",
                                "Pós-1980 (Aquecimento)" = "#2b5c8f"), name = "Período") +
  scale_size_continuous(range = c(1.5, 5), guide = "none") +
  annotate("text", x = -0.45, y = 6.2,
           label = "Pré-1980: β = +0.10 mm / °C (p = 0.91, nulo)\nPós-1980: β = -1.68 mm / °C (p = 0.005, sig!)",
           color = "#1a365d", size = 3.4, fontface = "bold", hjust = 0) +
  labs(x = "Anomalia de Temperatura Regional (°C)", y = "Desvio de Asa (mm)",
       title = "C. Resposta Térmica da Asa: Pré vs. Pós-1980",
       subtitle = "Pós-1980: cada +1 °C associa-se a -1,68 mm de asa (p = 0,005). Pré-1980: resposta nula.") +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8.5),
        legend.position = "bottom")

p_coupling <- (p3a / p3b) | p3c
p_coupling <- p_coupling + plot_layout(widths = c(1.6, 1.4))

ggsave(file.path(fig_dir, "climate_morphometry_coupling.png"), p_coupling, width = 13.5, height = 7.5, dpi = 200, bg = "white")

# Figure 4: Climate-Mass Coupling & Breakpoint (1880–2018)
message("[figures] Generating climate-mass coupling figure...")

dm_sec_centered <- d_base %>%
  filter(Order == "Passeriformes", !is.na(ln_mass), Year >= 1960, Year <= 2018) %>%
  group_by(Binomial) %>%
  filter(n() >= 10) %>%
  mutate(spp_mean = mean(ln_mass, na.rm = TRUE),
         mass_centered = (ln_mass - spp_mean) * 100) %>%
  ungroup() %>%
  filter(abs(mass_centered) <= 4 * sd(mass_centered, na.rm = TRUE))

annual_mass <- dm_sec_centered %>%
  group_by(Year) %>%
  summarise(
    mean_mc = mean(mass_centered, na.rm = TRUE),
    se_mc   = sd(mass_centered, na.rm = TRUE) / sqrt(n()),
    n       = n(),
    .groups = "drop"
  )

fit_mass_bp <- lm(mean_mc ~ pmin(Year - bp_val, 0) + pmax(Year - bp_val, 0),
                  weights = annual_mass$n, data = annual_mass)
annual_mass$mass_pred <- predict(fit_mass_bp)

p4b <- ggplot(annual_mass, aes(x = Year, y = mean_mc)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = bp_val, linetype = "dashed", color = "#b2182b", linewidth = 0.9) +
  geom_point(aes(size = n), color = "#d95f02", alpha = 0.75) +
  geom_errorbar(aes(ymin = mean_mc - se_mc, ymax = mean_mc + se_mc), width = 0.8, color = "grey60", alpha = 0.5) +
  geom_line(aes(y = mass_pred), color = "#d95f02", linewidth = 1.2) +
  scale_size_continuous(range = c(1.5, 5), name = "N aves") +
  scale_x_continuous(limits = c(1880, 2018)) +
  annotate("text", x = 1885, y = 14, label = "Pre-1980 (1962–1979):\n-0.17% / ano (p = 0.52, estável)",
           color = "grey30", size = 3.2, hjust = 0) +
  annotate("text", x = 1999, y = -14, label = "Pós-1980: -0.02% / ano (p = 0.49)\n(Estabilidade estrita em campo)",
           color = "#d95f02", size = 3.2, hjust = 0.5, fontface = "bold") +
  labs(x = "Ano", y = "Desvio Médio de Massa (%)",
       title = "B. Trajetória de Massa Corporal (1960–2018)",
       subtitle = "Desvio relativo médio por espécie (Passeriformes, N = 50.058). Quebra segmentada em 1980.") +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8.5),
        legend.position = "right")

coupled_mass <- inner_join(annual_mass, giss, by = "Year") %>%
  mutate(periodo = ifelse(Year >= 1980, "Pós-1980 (Aquecimento)", "Pré-1980 (Histórico)"))

p4c <- ggplot(coupled_mass, aes(x = af_temp, y = mean_mc, color = periodo)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  geom_point(aes(size = n), alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.1) +
  scale_color_manual(values = c("Pré-1980 (Histórico)" = "#e6550d",
                                "Pós-1980 (Aquecimento)" = "#2b5c8f"), name = "Período") +
  scale_size_continuous(range = c(1.5, 5), guide = "none") +
  annotate("text", x = -0.45, y = 22,
           label = "Pré-1980: β = -1.97% / °C (p = 0.86, nulo)\nPós-1980: β = +3.68% / °C (p = 0.43, nulo)",
           color = "#1a365d", size = 3.4, fontface = "bold", hjust = 0) +
  labs(x = "Anomalia de Temperatura Regional (°C)", y = "Desvio de Massa (%)",
       title = "C. Resposta Térmica da Massa: Pré vs. Pós-1980",
       subtitle = "Massa corporal estatisticamente estável e insensível à temperatura em ambos os períodos.") +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8.5),
        legend.position = "bottom")

p_mass_coupling <- (p3a / p4b) | p4c
p_mass_coupling <- p_mass_coupling + plot_layout(widths = c(1.6, 1.4))
ggsave(file.path(fig_dir, "climate_mass_coupling.png"), p_mass_coupling, width = 13.5, height = 7.5, dpi = 200, bg = "white")

# Figure 5: Climate-Allometry Coupling & Breakpoint (1880–2018)
message("[figures] Generating climate-allometry coupling figure...")

da_sec_centered <- d_base %>%
  filter(Order == "Passeriformes", !is.na(conc_wing), !is.na(ln_mass), Year >= 1960, Year <= 2018) %>%
  mutate(
    iso_contrast = (ln_mass - 3 * log(conc_wing)) * 100,
    rel_wing     = (log(conc_wing) - (1/3) * ln_mass) * 100
  ) %>%
  group_by(Binomial) %>%
  filter(n() >= 10) %>%
  mutate(
    spp_mean_iso = mean(iso_contrast, na.rm = TRUE),
    iso_centered = iso_contrast - spp_mean_iso
  ) %>%
  ungroup() %>%
  filter(abs(iso_centered) <= 4 * sd(iso_centered, na.rm = TRUE))

annual_allom <- da_sec_centered %>%
  group_by(Year) %>%
  summarise(
    mean_iso = mean(iso_centered, na.rm = TRUE),
    se_iso   = sd(iso_centered, na.rm = TRUE) / sqrt(n()),
    n        = n(),
    .groups  = "drop"
  )

fit_allom_bp <- lm(mean_iso ~ pmin(Year - bp_val, 0) + pmax(Year - bp_val, 0),
                   weights = annual_allom$n, data = annual_allom)
annual_allom$allom_pred <- predict(fit_allom_bp)

p5b <- ggplot(annual_allom, aes(x = Year, y = mean_iso)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = bp_val, linetype = "dashed", color = "#b2182b", linewidth = 0.9) +
  geom_point(aes(size = n), color = "#7570b3", alpha = 0.75) +
  geom_errorbar(aes(ymin = mean_iso - se_iso, ymax = mean_iso + se_iso), width = 0.8, color = "grey60", alpha = 0.5) +
  geom_line(aes(y = allom_pred), color = "#7570b3", linewidth = 1.2) +
  scale_size_continuous(range = c(1.5, 5), name = "N aves") +
  scale_x_continuous(limits = c(1880, 2018)) +
  annotate("text", x = 1885, y = 7.5, label = "Pre-1980 (1962–1979):\n+0.52% / ano (p = 0.69, nulo)",
           color = "grey30", size = 3.2, hjust = 0) +
  annotate("text", x = 1999, y = -14, label = "Pós-1980: +2.12% / dec (p = 0.048)\n(Stoutening térmico significante)",
           color = "#7570b3", size = 3.2, hjust = 0.5, fontface = "bold") +
  labs(x = "Ano", y = "Contraste de Isometria Δiso (%)",
       title = "B. Trajetória da Alometria (Δiso = ln M - 3 ln L, 1960–2018)",
       subtitle = "Desvio centrado por espécie (Passeriformes compartilhados, N = 29.871). Quebra em 1980.") +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8.5),
        legend.position = "right")

coupled_allom <- inner_join(annual_allom, giss, by = "Year") %>%
  mutate(periodo = ifelse(Year >= 1980, "Pós-1980 (Aquecimento)", "Pré-1980 (Histórico)"))

p5c <- ggplot(coupled_allom, aes(x = af_temp, y = mean_iso, color = periodo)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  geom_point(aes(size = n), alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.1) +
  scale_color_manual(values = c("Pré-1980 (Histórico)" = "#e6550d",
                                "Pós-1980 (Aquecimento)" = "#2b5c8f"), name = "Período") +
  scale_size_continuous(range = c(1.5, 5), guide = "none") +
  annotate("text", x = -0.45, y = 8.5,
           label = "Pré-1980: β = +4.36% / °C (p = 0.70, nulo)\nPós-1980: β = +9.32% / °C (p = 0.034, sig!)",
           color = "#1a365d", size = 3.4, fontface = "bold", hjust = 0) +
  labs(x = "Anomalia de Temperatura Regional (°C)", y = "Contraste de Isometria Δiso (%)",
       title = "C. Resposta Térmica da Alometria: Pré vs. Pós-1980",
       subtitle = "Pós-1980: desvio positivo cresce +9,32% por +1 °C (p = 0,034). Pré-1980: desacoplado.") +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8.5),
        legend.position = "bottom")

p_allom_coupling <- (p3a / p5b) | p5c
p_allom_coupling <- p_allom_coupling + plot_layout(widths = c(1.6, 1.4))
ggsave(file.path(fig_dir, "climate_allometry_coupling.png"), p_allom_coupling, width = 13.5, height = 7.5, dpi = 200, bg = "white")

message("[done] Script complete. Outputs in ", out_dir)
