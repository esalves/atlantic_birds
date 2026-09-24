// phylo_lmm_marginal.stan
// ---------------------------------------------------------------------------
// Gaussian phylogenetic mixed model for the Atlantic-bird record-level data,
// written to replace the brms fits (hours per tree) with a program that samples
// in minutes. Adapted from the marginalized spectral models of the
// Mammalian_decomposition project (S. Ortega; stan_production/), whose trick is
// to integrate Gaussian latent effects out of the likelihood analytically.
//
// That project's exact marginalization needs one observation per species. Our
// data are individual captures with crossed contributor / site / individual
// effects, so the species effects cannot be integrated out. What CAN be
// integrated out exactly is the individual (ring) random intercept, which is
// the largest block of latent parameters (thousands of rings) and the main
// cause of slow, funnel-shaped brms posteriors. Records of one individual have
// covariance D_i + tau^2 11', whose inverse and determinant have closed forms
// (Sherman-Morrison), with or without a heteroscedastic residual D_i. The
// individual effects are therefore never sampled.
//
// Everything else follows the brms/glmmTMB model:
//   y = X beta + a[spp] + b[spp] * t + sum_j c_j[group_j] + sum_k d_k[sgroup_k] * x_k
//       + u[ind] + e,            e ~ N(0, exp(W gamma + sum_k f_k[dgroup_k])^2)
// Species intercepts a and slopes b each split into a phylogenetic and an iid
// part through a total SD and a heritability-like share h2 (Santiago's h2
// reparameterization), with the phylogenetic part drawn in the eigenbasis of A
// (LA = U diag(sqrt(lambda)), so LA LA' = A):
//   a = sd_a * (sqrt(h2_a) * LA z_ap + sqrt(1 - h2_a) * z_ai)
//   b = sd_b * (sqrt(h2_b) * LA z_bp + sqrt(1 - h2_b) * z_bi)
// Setting phylo_slope = 1 gives the species slopes a Brownian component; with
// phylo_cor = 1 the phylogenetic intercept and slope are correlated
// (covariance Sigma (x) A), which glmmTMB's propto cannot fit.
//
// y is standardized in R before fitting; R back-transforms the draws.
// ---------------------------------------------------------------------------
data {
  int<lower=1> N;
  vector[N] y;
  int<lower=1> K;
  matrix[N, K] X;                       // fixed effects, first column intercept

  // species
  int<lower=1> S;
  array[N] int<lower=1, upper=S> spp;
  int<lower=0, upper=1> has_slope;      // species random slope on t
  vector[N] t;
  matrix[S, S] LA;                      // U diag(sqrt(lambda)) of the correlation matrix A
  int<lower=0, upper=1> phylo_int;
  int<lower=0, upper=1> phylo_slope;
  int<lower=0, upper=1> phylo_cor;

  // random-intercept groups (contributor, site, locality-year ...)
  int<lower=0> J;
  array[J] int<lower=1> Gj;
  array[N, J] int<lower=1> gidx;
  // random-slope groups (e.g. (0 + scaled_yr | src) for M4)
  int<lower=0> Js;
  array[Js] int<lower=1> Gs;
  array[N, Js] int<lower=1> sidx;
  matrix[N, Js] xs;

  // individual effect, marginalized: CSR incidence matrix (n_ind x N) of ones
  int<lower=0, upper=1> has_ind;
  int<lower=0> n_ind;
  array[has_ind ? N : 0] int<lower=1> csr_v;
  array[has_ind ? n_ind + 1 : 0] int<lower=1> csr_u;

  // residual (log-sigma) model
  int<lower=1> P;
  matrix[N, P] W;                       // first column intercept
  int<lower=0> Jd;                      // random intercepts in log-sigma (variance ladder)
  array[Jd] int<lower=1> Gd;
  array[N, Jd] int<lower=1> didx;

  real<lower=0> prior_beta_sd;
  real<lower=0> prior_sd_scale;
}

transformed data {
  vector[has_ind ? N : 0] csr_w = rep_vector(1.0, has_ind ? N : 0);
  int n_gj = sum(Gj);
  int n_gs = sum(Gs);
  int n_gd = sum(Gd);
}

parameters {
  vector[K] beta;
  vector[P] gamma;

  real<lower=0> sd_a;
  array[phylo_int] real<lower=0, upper=1> h2_a;
  vector[S] z_ai;
  vector[phylo_int ? S : 0] z_ap;

  array[has_slope] real<lower=0> sd_b;
  array[has_slope * phylo_slope] real<lower=0, upper=1> h2_b;
  vector[has_slope ? S : 0] z_bi;
  vector[has_slope * phylo_slope ? S : 0] z_bp;
  cholesky_factor_corr[phylo_cor ? 2 : 1] L_R;   // 1 x 1 (no parameter) unless phylo_cor = 1

  vector<lower=0>[J] sd_g;
  vector[n_gj] z_g;
  vector<lower=0>[Js] sd_s;
  vector[n_gs] z_s;

  vector<lower=0>[Jd] sd_d;
  vector[n_gd] z_d;

  array[has_ind] real<lower=0> tau;
}

transformed parameters {
  vector[S] a;
  vector[S] b = rep_vector(0, S);
  {
    if (phylo_int == 1)
      a = sd_a * (sqrt(h2_a[1]) * (LA * z_ap) + sqrt(1 - h2_a[1]) * z_ai);
    else
      a = sd_a * z_ai;
    if (has_slope == 1) {
      if (phylo_slope == 1) {
        vector[S] zb = z_bp;
        if (phylo_cor == 1) zb = L_R[2, 1] * z_ap + L_R[2, 2] * z_bp;
        b = sd_b[1] * (sqrt(h2_b[1]) * (LA * zb) + sqrt(1 - h2_b[1]) * z_bi);
      } else {
        b = sd_b[1] * z_bi;
      }
    }
  }
}

model {
  vector[N] eta = X * beta + a[spp];
  vector[N] log_sd = W * gamma;
  if (has_slope == 1) eta += b[spp] .* t;
  {
    int pos = 1;
    for (j in 1:J) {
      vector[Gj[j]] c = sd_g[j] * segment(z_g, pos, Gj[j]);
      eta += c[gidx[, j]];
      pos += Gj[j];
    }
    pos = 1;
    for (k in 1:Js) {
      vector[Gs[k]] d = sd_s[k] * segment(z_s, pos, Gs[k]);
      eta += d[sidx[, k]] .* xs[, k];
      pos += Gs[k];
    }
    pos = 1;
    for (k in 1:Jd) {
      vector[Gd[k]] e = sd_d[k] * segment(z_d, pos, Gd[k]);
      log_sd += e[didx[, k]];
      pos += Gd[k];
    }
  }

  // priors (y is standardized)
  beta ~ normal(0, prior_beta_sd);
  gamma ~ normal(0, 1);
  sd_a ~ student_t(3, 0, prior_sd_scale);
  sd_b ~ student_t(3, 0, prior_sd_scale);
  sd_g ~ student_t(3, 0, prior_sd_scale);
  sd_s ~ student_t(3, 0, prior_sd_scale);
  sd_d ~ student_t(3, 0, 1);
  z_d ~ std_normal();
  tau ~ student_t(3, 0, prior_sd_scale);
  h2_a ~ beta(1, 1);
  h2_b ~ beta(1, 1);
  L_R ~ lkj_corr_cholesky(2);
  z_ai ~ std_normal();
  z_ap ~ std_normal();
  z_bi ~ std_normal();
  z_bp ~ std_normal();
  z_g ~ std_normal();
  z_s ~ std_normal();

  if (has_ind == 1) {
    // Sum over individuals of the exact log density of N(0, D_i + tau^2 11'):
    //   -0.5 [ sum log d + log(1 + tau^2 s0) + sum r^2/d - tau^2 s1^2 / (1 + tau^2 s0) ]
    // with s0 = sum_i 1/d, s1 = sum_i r/d.
    vector[N] d = exp(2 * log_sd);
    vector[N] r = y - eta;
    vector[N] rd = r ./ d;
    vector[n_ind] s0 = csr_matrix_times_vector(n_ind, N, csr_w, csr_v, csr_u, inv(d));
    vector[n_ind] s1 = csr_matrix_times_vector(n_ind, N, csr_w, csr_v, csr_u, rd);
    real t2 = square(tau[1]);
    vector[n_ind] q = 1 + t2 * s0;
    target += -0.5 * (N * log(2 * pi()) + 2 * sum(log_sd) + sum(log(q))
                      + dot_product(r, rd) - t2 * sum(square(s1) ./ q));
  } else {
    y ~ normal(eta, exp(log_sd));
  }
}

generated quantities {
  real cor_phylo_ab = 0;
  if (phylo_cor == 1) cor_phylo_ab = L_R[2, 1];
}
