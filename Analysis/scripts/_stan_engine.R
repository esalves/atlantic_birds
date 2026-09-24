# _stan_engine.R
# ---------------------------------------------------------------------------
# Bayesian counterpart of _phylo_engine.R: fits the same lme4-style model
# formulas with the custom Stan program stan/phylo_lmm_marginal.stan, one tree
# at a time, and pools the posterior draws across trees (the mixture posterior;
# Nakagawa & de Villemereuil 2019).
#
# WHY: brms took ~55 min per model per tree on Totoro, mostly because it samples
#      every ring (individual) intercept. The Stan program integrates those out
#      exactly and draws the phylogenetic effects in the eigenbasis of A, after
#      the marginalized spectral models of S. Ortega (Mammalian_decomposition,
#      stan_production/). See the header of the .stan file for the math.
#
# FORMULAS understood (the vocabulary of atlantic_parallel_controlled.R and
# referee_reruns.R):
#   (1 + x || spp)  species intercept + slope on x; phylogenetic parts are added
#                   by phylo_int / phylo_slope / phylo_cor
#   (1 | ind)       individual intercept -> marginalized (ind or ring; set by ind_group)
#   (1 | g)         any other random intercept
#   (1 + x || g)    random intercept + independent random slope on x
# Correlated non-species terms (1 + x | g) are not supported.
# dispformula becomes the log-sigma model (Phase 3 variance ladder).
#
# BACKEND: cmdstanr when installed (Totoro), rstan otherwise (laptop).
# ---------------------------------------------------------------------------
suppressMessages({ library(dplyr) })
`%||%` <- function(a, b) if (is.null(a)) b else a

STAN_FILE <- file.path(if (exists("ANALYSIS_DIR")) ANALYSIS_DIR else ".", "scripts", "stan", "phylo_lmm_marginal.stan")

.stan_backend <- function() if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"

#' Compile once per session.
stan_model_cached <- local({
  mod <- NULL
  function(file = STAN_FILE) {
    if (!is.null(mod)) return(mod)
    mod <<- if (.stan_backend() == "cmdstanr") cmdstanr::cmdstan_model(file)
            else rstan::stan_model(file, auto_write = TRUE)
    mod
  }
})

#' LA = U diag(sqrt(lambda)) so that LA %*% t(LA) == A (spectral square root).
spectral_root <- function(A) {
  e <- eigen(A, symmetric = TRUE)
  L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)))
  dimnames(L) <- list(rownames(A), NULL); L
}

#' Translate formula + data (+ A) into the Stan data list.
make_stan_data <- function(formula, data, A, dispformula = ~1, species_group = "spp",
                           species_col = "species_name", ind_group = c("ind", "ring"),
                           phylo_int = TRUE, phylo_slope = FALSE, phylo_cor = FALSE,
                           prior_beta_sd = 5, prior_sd_scale = 1) {
  if (phylo_cor && !(phylo_int && phylo_slope)) stop("phylo_cor needs phylo_int and phylo_slope")
  bars  <- lme4::findbars(lme4::expandDoubleVerts(formula))
  fixed <- lme4::nobars(formula)
  vars  <- unique(c(all.vars(formula), all.vars(dispformula), species_col))
  data  <- as.data.frame(data)[, intersect(vars, names(data)), drop = FALSE]
  data  <- data[stats::complete.cases(data), , drop = FALSE]
  sp <- as.character(data[[species_col]])
  data <- data[sp %in% rownames(A), , drop = FALSE]
  sp_f <- factor(as.character(data[[species_col]]), levels = rownames(A))
  if (any(table(sp_f) == 0)) {                    # keep A to the species present
    keep <- names(which(table(sp_f) > 0)); A <- A[keep, keep]; sp_f <- factor(as.character(sp_f), levels = keep)
  }

  y_raw <- model.response(model.frame(fixed, data))
  y_m <- mean(y_raw); y_s <- sd(y_raw)
  X <- model.matrix(fixed, data)
  dbars <- lme4::findbars(dispformula)
  W <- model.matrix(lme4::nobars(dispformula) %||% ~1, data)
  D_idx <- list()
  for (b in dbars) {
    if (deparse(b[[2]]) != "1") stop("only random intercepts are supported in dispformula: ", deparse(b))
    g <- deparse(b[[3]])
    D_idx[[g]] <- as.integer(factor(interaction(lapply(strsplit(g, ":", fixed = TRUE)[[1]], function(v) data[[v]]), drop = TRUE)))
  }
  stopifnot(colnames(X)[1] == "(Intercept)", colnames(W)[1] == "(Intercept)")

  slope_var <- NULL; J_idx <- list(); S_idx <- list(); S_x <- list(); ind <- NULL
  for (b in bars) {
    g   <- deparse(b[[3]])
    lhs <- attr(terms(as.formula(paste("~", deparse(b[[2]])))), "term.labels")
    has_int <- attr(terms(as.formula(paste("~", deparse(b[[2]])))), "intercept") == 1
    if (length(lhs) > 1) stop("correlated random slopes not supported: ", deparse(b))
    gf <- factor(interaction(lapply(strsplit(g, ":", fixed = TRUE)[[1]], function(v) data[[v]]), drop = TRUE))
    if (g == species_group) {
      if (length(lhs)) slope_var <- lhs
      next
    }
    if (g %in% ind_group && has_int && !length(lhs) && is.null(ind)) { ind <- gf; next }
    if (has_int && !length(lhs)) J_idx[[g]] <- as.integer(gf)
    else if (!has_int && length(lhs) == 1) { S_idx[[paste0(g, ":", lhs)]] <- as.integer(gf); S_x[[paste0(g, ":", lhs)]] <- data[[lhs]] }
    else stop("unsupported term: ", deparse(b))
  }

  N <- nrow(data)
  sd_list <- list(
    N = N, y = (y_raw - y_m) / y_s, K = ncol(X), X = X,
    S = nlevels(sp_f), spp = as.integer(sp_f),
    has_slope = as.integer(!is.null(slope_var)),
    t = if (is.null(slope_var)) rep(0, N) else data[[slope_var]],
    LA = spectral_root(A),
    phylo_int = as.integer(phylo_int), phylo_slope = as.integer(phylo_slope && !is.null(slope_var)),
    phylo_cor = as.integer(phylo_cor),
    J = length(J_idx), Gj = as.array(vapply(J_idx, max, 1L)),
    gidx = if (length(J_idx)) do.call(cbind, J_idx) else matrix(1L, N, 0),
    Js = length(S_idx), Gs = as.array(vapply(S_idx, max, 1L)),
    sidx = if (length(S_idx)) do.call(cbind, S_idx) else matrix(1L, N, 0),
    xs = if (length(S_x)) do.call(cbind, S_x) else matrix(0, N, 0),
    has_ind = as.integer(!is.null(ind)),
    P = ncol(W), W = W,
    Jd = length(D_idx), Gd = as.array(vapply(D_idx, max, 1L)),
    didx = if (length(D_idx)) do.call(cbind, D_idx) else matrix(1L, N, 0),
    prior_beta_sd = prior_beta_sd, prior_sd_scale = prior_sd_scale
  )
  if (!is.null(ind)) {                               # CSR incidence, rows = individuals
    ii <- as.integer(droplevels(ind)); ord <- order(ii)
    sd_list$n_ind <- max(ii); sd_list$csr_v <- ord
    sd_list$csr_u <- c(1L, cumsum(tabulate(ii, max(ii))) + 1L)
  } else { sd_list$n_ind <- 0L; sd_list$csr_v <- integer(0); sd_list$csr_u <- integer(0) }
  if (!length(J_idx)) sd_list$Gj <- integer(0)
  if (!length(S_idx)) sd_list$Gs <- integer(0)
  if (!length(D_idx)) sd_list$Gd <- integer(0)

  attr(sd_list, "meta") <- list(y_mean = y_m, y_sd = y_s, X_names = colnames(X), W_names = colnames(W),
                                J_names = names(J_idx), Js_names = names(S_idx), Jd_names = names(D_idx), species = levels(sp_f),
                                slope_var = slope_var, has_ind = !is.null(ind), N = N)
  sd_list
}

#' Back-transform draws to the response scale and label them.
backtransform_draws <- function(dr, meta) {
  s <- meta$y_sd; out <- list()
  B <- as.matrix(dr[, grep("^beta\\[", names(dr)), drop = FALSE]) * s
  B[, 1] <- B[, 1] + meta$y_mean; colnames(B) <- paste0("b_", meta$X_names)
  G <- as.matrix(dr[, grep("^gamma\\[", names(dr)), drop = FALSE]); G[, 1] <- G[, 1] + log(s)
  colnames(G) <- paste0("sigma_", meta$W_names)
  sds <- cbind(sd_spp_int = dr[["sd_a"]] * s)
  if ("sd_b[1]" %in% names(dr)) sds <- cbind(sds, sd_spp_slope = dr[["sd_b[1]"]] * s)
  if ("h2_a[1]" %in% names(dr)) sds <- cbind(sds, h2_int = dr[["h2_a[1]"]])
  if ("h2_b[1]" %in% names(dr)) sds <- cbind(sds, h2_slope = dr[["h2_b[1]"]])
  add <- function(m, nm, v) { m <- cbind(m, v); colnames(m)[ncol(m)] <- nm; m }
  for (j in seq_along(meta$J_names)) sds <- add(sds, paste0("sd_", meta$J_names[j]), dr[[sprintf("sd_g[%d]", j)]] * s)
  for (k in seq_along(meta$Js_names)) sds <- add(sds, paste0("sd_", meta$Js_names[k]), dr[[sprintf("sd_s[%d]", k)]] * s)
  for (k in seq_along(meta$Jd_names)) sds <- add(sds, paste0("sdlogsigma_", meta$Jd_names[k]), dr[[sprintf("sd_d[%d]", k)]])
  if (meta$has_ind) sds <- add(sds, "sd_ind", dr[["tau[1]"]] * s)
  sds <- cbind(sds, cor_phylo_int_slope = dr[["cor_phylo_ab"]])
  bsp <- as.matrix(dr[, grep("^b\\[", names(dr)), drop = FALSE]) * s
  colnames(bsp) <- paste0("slope_dev_", meta$species)
  cbind(B, G, sds, bsp)
}

#' Fit one tree. Returns list(draws (matrix, response scale), diag, secs, meta).
fit_stan_tree <- function(formula, data, A, dispformula = ~1, chains = 4, iter_warmup = 1000,
                          iter_sampling = 500, adapt_delta = 0.95, max_treedepth = 12, seed = 20240101,
                          parallel_chains = chains, ...) {
  sdat <- make_stan_data(formula, data, A, dispformula = dispformula, ...)
  meta <- attr(sdat, "meta"); mod <- stan_model_cached()
  keep <- c("beta", "gamma", "sd_a", "h2_a", "sd_b", "h2_b", "sd_g", "sd_s", "sd_d", "tau", "b", "cor_phylo_ab")
  t0 <- Sys.time()
  if (.stan_backend() == "cmdstanr") {
    fit <- mod$sample(data = sdat, chains = chains, parallel_chains = parallel_chains,
                      iter_warmup = iter_warmup, iter_sampling = iter_sampling, seed = seed,
                      adapt_delta = adapt_delta, max_treedepth = max_treedepth, refresh = 0, show_messages = FALSE)
    dr <- posterior::as_draws_df(fit$draws(variables = intersect(keep, fit$metadata()$stan_variables)))
    sm <- fit$diagnostic_summary(quiet = TRUE)
    diag <- list(divergent = sum(sm$num_divergent), treedepth_hits = sum(sm$num_max_treedepth),
                 ebfmi_min = min(sm$ebfmi))
  } else {
    fit <- rstan::sampling(mod, data = sdat, chains = chains, cores = parallel_chains,
                           warmup = iter_warmup, iter = iter_warmup + iter_sampling, seed = seed,
                           control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth), refresh = 0)
    dr <- posterior::as_draws_df(as.array(fit, pars = intersect(keep, fit@model_pars)))
    sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    diag <- list(divergent = sum(sapply(sp, function(x) sum(x[, "divergent__"]))),
                 treedepth_hits = sum(sapply(sp, function(x) sum(x[, "treedepth__"] >= max_treedepth))),
                 ebfmi_min = NA_real_)
  }
  secs <- as.numeric(Sys.time() - t0, units = "secs")
  summ <- posterior::summarise_draws(posterior::subset_draws(dr, variable = c("beta", "sd_a")), "rhat", "ess_bulk")
  diag$rhat_max <- max(summ$rhat, na.rm = TRUE); diag$ess_bulk_min <- min(summ$ess_bulk, na.rm = TRUE)
  draws <- backtransform_draws(as.data.frame(dr), meta)
  list(draws = draws, chain = dr$.chain, diag = diag, secs = secs, meta = meta)
}

#' Pool per-tree draw matrices (mixture posterior) into a summary table.
pool_stan_draws <- function(per_tree, pars = NULL) {
  D <- do.call(rbind, lapply(per_tree, `[[`, "draws"))
  if (!is.null(pars)) D <- D[, intersect(pars, colnames(D)), drop = FALSE]
  data.frame(par = colnames(D), estimate = colMeans(D), sd = apply(D, 2, sd),
             lower = apply(D, 2, quantile, 0.025), upper = apply(D, 2, quantile, 0.975),
             p_neg = colMeans(D < 0), n_draws = nrow(D), n_trees = length(per_tree), row.names = NULL)
}

#' Loop over trees (serial; on Totoro the per-tree driver parallelizes instead).
run_phylo_trees_stan <- function(formula, data, A_list, trees = seq_along(A_list), verbose = TRUE, ...) {
  res <- lapply(trees, function(i) {
    r <- fit_stan_tree(formula, data, A_list[[i]], ...); r$tree <- i
    if (verbose) message(sprintf("[stan] tree %d: %.0f s, rhat_max %.3f, ess_min %.0f, divergent %d",
                                 i, r$secs, r$diag$rhat_max, r$diag$ess_bulk_min, r$diag$divergent))
    r
  })
  list(pooled = pool_stan_draws(res),
       diag = do.call(rbind, lapply(res, function(r) data.frame(tree = r$tree, secs = r$secs, r$diag))),
       per_tree = res)
}
