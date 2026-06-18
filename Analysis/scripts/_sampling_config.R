# _sampling_config.R
# ---------------------------------------------------------------------------
# Shared MCMC settings for the phylogenetic brms models, sourced by
# atlantic_parallel.R, atlantic_parallel_bill.R, atlantic_parallel_mass.R,
# atlantic_diet_interaction.R and atlantic_sensitivity_thresholds.R.
#
# WHY these statistical settings: the original single-chain, iter=1500/warmup=1000
# left the POPULATION INTERCEPT poorly mixed (R-hat up to ~1.07, bulk-ESS ~11) —
# weak identifiability between the fixed intercept and the phylogenetic group
# intercept when phylogenetic signal ~0.95. The fix is more/longer/multiple chains
# with tighter adaptation: chains 4, iter 4000 (warmup 2000), adapt_delta 0.99,
# max_treedepth 15. These STATISTICAL settings are FIXED across hosts.
#
# HOST-AWARE CONCURRENCY (added 2026-06): the file auto-detects whether it is on
# Totoro (big shared server) vs. a laptop, and sets only the *concurrency* and
# *backend* accordingly, so the SAME file is safe on both:
#   * Totoro  -> cmdstanr, 4 chains in parallel (cores=4), 46 tree-workers
#                (46*4 = 184 cores/job, so two model jobs fit within 384 cores).
#   * laptop  -> rstan, chains sequential (cores=1), ONE tree-worker by default.
# The laptop hazard is that each worker holds a full brms fit in memory at once;
# running one-per-core OOMs a laptop. Off-server we default to a SINGLE worker.
# Bump it deliberately with: ATLANTIC_WORKERS=2 Rscript <script>.R
# Force the server profile anywhere with: ATLANTIC_PROFILE=server
# ---------------------------------------------------------------------------

.node    <- tolower(Sys.info()[["nodename"]])
.profile <- Sys.getenv("ATLANTIC_PROFILE", "")
.on_server <- if (nzchar(.profile)) identical(.profile, "server") else grepl("totoro", .node)

SAMPLING <- list(
  chains        = 4L,
  iter          = 4000L,
  warmup        = 2000L,
  adapt_delta   = 0.99,
  max_treedepth = 15L,
  # cmdstanr on the server (faster); rstan elsewhere (the portable default).
  backend       = if (.on_server && requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan",
  # chains in parallel only on the server; sequential on a laptop.
  cores         = if (.on_server) 4L else 1L
)

# Across-tree future_lapply concurrency.
SAMPLING$workers <- if (.on_server) {
  46L
} else {
  n <- suppressWarnings(as.integer(Sys.getenv("ATLANTIC_WORKERS", "1")))
  if (is.na(n) || n < 1L) 1L else n   # default 1 worker off-server = OOM-safe
}

# Convenience: per-fit control list, so scripts can pass `control = SAMPLING_CONTROL`.
SAMPLING_CONTROL <- list(adapt_delta = SAMPLING$adapt_delta,
                         max_treedepth = SAMPLING$max_treedepth)

message(sprintf("[sampling] host=%s chains=%d iter=%d warmup=%d adapt_delta=%.2f backend=%s cores=%d workers=%d",
                if (.on_server) "server" else "local",
                SAMPLING$chains, SAMPLING$iter, SAMPLING$warmup,
                SAMPLING$adapt_delta, SAMPLING$backend, SAMPLING$cores, SAMPLING$workers))
