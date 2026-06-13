# _sampling_config.R
# ---------------------------------------------------------------------------
# Shared MCMC settings for the phylogenetic brms models, sourced by
# atlantic_parallel.R, atlantic_parallel_bill.R, atlantic_parallel_mass.R,
# atlantic_diet_interaction.R and atlantic_sensitivity_thresholds.R.
#
# WHY: the original single-chain, iter=1500/warmup=1000 settings left the
# POPULATION INTERCEPT poorly mixed (R-hat up to ~1.07, bulk-ESS as low as ~11).
# That is the classic weak identifiability between the fixed intercept and the
# phylogenetic group intercept when phylogenetic signal is ~0.95 — almost all of
# the among-species variance is phylogenetic, so the two intercepts trade off.
# The fix is more, longer, multiple chains with tighter adaptation:
#   * chains 1 -> 4        proper (multi-chain) R-hat + 4x the draws
#   * warmup 1000 -> 2000  longer adaptation for the funnel-shaped geometry
#   * iter   1500 -> 4000  (=> 2000 post-warmup draws/chain, 8000 total)
#   * adapt_delta 0.95 -> 0.99, max_treedepth 12 -> 15  (fewer divergences)
# These are heavier but intended for a many-core server (Totoro). EDIT THIS FILE
# to retune for the run host — every model script picks the values up from here.
# ---------------------------------------------------------------------------

SAMPLING <- list(
  chains        = 4L,
  iter          = 4000L,
  warmup        = 2000L,
  adapt_delta   = 0.99,
  max_treedepth = 15L,
  # Stan backend. "rstan" is the safe default (matches the current toolchain).
  # On Totoro, if a cmdstan toolchain is installed, set this to "cmdstanr" for a
  # substantial speed-up (and optionally enable within-chain threading below).
  backend       = "rstan",
  # cores per fit. Chains run SEQUENTIALLY within each fit (cores = 1) so the
  # outer per-tree future parallelism does not nest into brms's per-chain
  # parallelism (which can oversubscribe or hang). With workers = physical cores
  # the many tree fits already saturate the machine. Raise only if you reduce
  # `workers` correspondingly so workers * cores <= physical cores.
  cores         = 1L
)

# Concurrency for the across-tree future_lapply: one worker per physical core.
SAMPLING$workers <- local({
  n <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
  if (is.na(n) || n < 1L) 4L else as.integer(n)
})

# Convenience: the per-fit control list, so scripts can pass `control = SAMPLING_CONTROL`.
SAMPLING_CONTROL <- list(adapt_delta = SAMPLING$adapt_delta,
                         max_treedepth = SAMPLING$max_treedepth)

message(sprintf("[sampling] chains=%d iter=%d warmup=%d adapt_delta=%.2f backend=%s workers=%d",
                SAMPLING$chains, SAMPLING$iter, SAMPLING$warmup,
                SAMPLING$adapt_delta, SAMPLING$backend, SAMPLING$workers))
