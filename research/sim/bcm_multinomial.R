## sim/bcm_multinomial.R
##
## Empirical-likelihood test of beta in the multinomial choice model
##
##   U_{ij} = beta_j^T x_i + eps_{ij},  j = 1, ..., J-1,   U_{i0} = 0,
##   y_i = argmax_j U_{ij}.
##
## Strategy: each pair (j, 0) gives a standard BCM with coefficient
## gamma_j = beta_j - beta_0 and residual eps_{i0} - eps_{ij}.  We
## test H0: gamma = (gamma_1^0, ..., gamma_{J-1}^0) by stacking the
## J-1 pairwise BJ-estimating functions and running Owen-EL with one
## scale coordinate dropped.
##
## Under H0 and regularity (the BCM Wilks conditions of Section 6
## extended pairwise), the statistic is asymptotically
## chi^2_{p(J-1) - 1}.

source("sim/bcm_sim.R", chdir = FALSE)

## ---------------------------------------------------------
## 1. Pairwise BJ estimating function for one pair (j, 0)
##    Returns g matrix of dim (#subjects in pair) x p.
## ---------------------------------------------------------
pairwise_bj <- function(X_pair, delta_pair, gamma) {
  ## delta = 1 if subject chose alternative j; 0 if chose alt 0.
  bj <- bj_residuals(X_pair, delta_pair, gamma)
  bj$g
}

## ---------------------------------------------------------
## 2. Stack pairwise estimating functions across J-1 alternatives
##    Each row i of the stack has length (J-1)*p; padded with 0
##    in positions corresponding to pairs the subject doesn't belong to.
## ---------------------------------------------------------
stack_multinomial_g <- function(X, y, gammas) {
  ## y in {0, 1, ..., J-1};  gammas[[j]] is the BCM-coefficient for pair (j, 0).
  ## In BCM convention, delta = 1 means "chose the base (alt 0)", so
  ## delta_{i,j} = 1{y_i = 0} for pair (j, 0).  With this choice,
  ## delta = 1 iff the residual difference eta = eps_0 - eps_j satisfies
  ## eta > gamma_j^T x, which is the BCM form
  ## delta = 1{ -eta <= -gamma^T x } once we set tilde_eps = -eta.
  n <- nrow(X); p <- ncol(X); Jm <- length(gammas)
  G <- matrix(0, n, Jm * p)
  for (j in seq_len(Jm)) {
    idx <- which(y == j | y == 0L)
    if (length(idx) < 5L) next
    delta_pair <- as.integer(y[idx] == 0L)
    g_pair     <- pairwise_bj(X[idx, , drop = FALSE], delta_pair, gammas[[j]])
    cols <- ((j - 1L) * p + 1L):(j * p)
    G[idx, cols] <- g_pair
  }
  G
}

## ---------------------------------------------------------
## 3. Multinomial BCM-EL test statistic
## ---------------------------------------------------------
mn_bcm_lr <- function(X, y, gammas) {
  G <- stack_multinomial_g(X, y, gammas)
  ## Scale identification: each pair (j, 0) has its own BCM scale degeneracy.
  ## We drop the first coordinate of EACH pair's slot --- one drop per pair.
  ## After dropping, the system has rank (p - 1) * (J - 1).
  p   <- ncol(X)
  Jm  <- length(gammas)
  drop_cols <- ((seq_len(Jm) - 1L) * p) + 1L
  keep_cols <- setdiff(seq_len(Jm * p), drop_cols)
  G <- G[, keep_cols, drop = FALSE]
  keep <- rowSums(abs(G)) > 0
  G <- G[keep, , drop = FALSE]
  res <- tryCatch(
    emplik::el.test(G, mu = rep(0, ncol(G))),
    error = function(e) list(`-2LLR` = NA_real_)
  )
  list(stat = res[["-2LLR"]], df = ncol(G), n_eff = nrow(G))
}

## ---------------------------------------------------------
## 4. Data-generating process for multinomial (J alternatives)
##    eps_{ij} iid Gumbel (gives MNL model) or iid Logistic, etc.
## ---------------------------------------------------------
gen_multinomial <- function(n, betas, eps_dist = "gumbel") {
  ## betas: p x (J-1) matrix; column j gives beta_j
  p <- nrow(betas); Jm <- ncol(betas); J <- Jm + 1L
  X <- matrix(runif(n * p, -1, 1), n, p)
  U <- matrix(0, n, J)
  ## U[, 1] = 0 (base, j = 0)
  for (j in seq_len(Jm)) {
    U[, j + 1L] <- X %*% betas[, j]
  }
  ## add iid errors
  eps <- switch(eps_dist,
    gumbel   = -log(-log(matrix(runif(n * J), n, J))),  # Gumbel(0,1)
    logistic = matrix(rlogis(n * J), n, J),
    normal   = matrix(rnorm (n * J), n, J)
  )
  U <- U + eps
  y <- max.col(U) - 1L                # y in {0, 1, ..., J-1}
  list(X = X, y = y)
}

## ---------------------------------------------------------
## 5. One MC cell
## ---------------------------------------------------------
run_mn_cell <- function(n, betas, dgp_args, reps, seed) {
  set.seed(seed)
  Jm <- ncol(betas)
  gammas <- lapply(seq_len(Jm), function(j) betas[, j])  # since beta_0 = 0
  stats  <- numeric(reps)
  times  <- numeric(reps)
  df_use <- NA_integer_
  for (r in seq_len(reps)) {
    d  <- do.call(gen_multinomial, c(list(n = n, betas = betas), dgp_args))
    t0 <- proc.time()
    out <- mn_bcm_lr(d$X, d$y, gammas)
    times[r] <- (proc.time() - t0)[3L]
    stats[r] <- out$stat
    df_use   <- out$df
  }
  list(stats = stats, times = times, df = df_use)
}

summarise_mn <- function(cell) {
  s    <- cell$stats
  df   <- cell$df
  good <- !is.na(s) & is.finite(s) & s >= 0
  s    <- s[good]
  list(
    df          = df,
    n_ok        = sum(good),
    type1_10    = mean(s > qchisq(0.90, df)),
    type1_05    = mean(s > qchisq(0.95, df)),
    type1_01    = mean(s > qchisq(0.99, df)),
    mean_stat   = mean(s),
    median_stat = median(s),
    ks_p        = suppressWarnings(ks.test(s, "pchisq", df = df))$p.value,
    median_time = median(cell$times[good])
  )
}

## ---------------------------------------------------------
## 6. Main driver
## ---------------------------------------------------------
main_multinomial <- function(reps     = 2000L,
                             out_path = "sim/results/multinomial_results.rds") {
  cat(sprintf("\nMultinomial BCM-EL: reps=%d per cell\n", reps))
  cat(strrep("=", 70), "\n", sep = "")

  ## p=2 covariates, J=3 alternatives ⇒ df = p*(J-1) - 1 = 3
  betas_p2_J3 <- cbind(c(1.0, 0.5), c(-0.5, 1.0))     # 2 x 2 matrix
  ## p=2 covariates, J=4 alternatives ⇒ df = p*(J-1) - 1 = 5
  betas_p2_J4 <- cbind(c(1.0, 0.5), c(-0.5, 1.0), c(0.7, -0.3))  # 2 x 3
  ## p=3 covariates, J=3 alternatives ⇒ df = p*(J-1) - 1 = 5
  betas_p3_J3 <- cbind(c(1.0, 0.5, -0.3), c(-0.5, 1.0, 0.2))

  cfgs <- list(
    list(name = "gumbel_p2_J3_n500",  n = 500L,  betas = betas_p2_J3, dgp_args = list(eps_dist = "gumbel")),
    list(name = "gumbel_p2_J3_n1000", n = 1000L, betas = betas_p2_J3, dgp_args = list(eps_dist = "gumbel")),
    list(name = "gumbel_p2_J3_n2000", n = 2000L, betas = betas_p2_J3, dgp_args = list(eps_dist = "gumbel")),
    list(name = "normal_p2_J3_n1000", n = 1000L, betas = betas_p2_J3, dgp_args = list(eps_dist = "normal")),
    list(name = "gumbel_p2_J4_n1000", n = 1000L, betas = betas_p2_J4, dgp_args = list(eps_dist = "gumbel")),
    list(name = "gumbel_p3_J3_n1000", n = 1000L, betas = betas_p3_J3, dgp_args = list(eps_dist = "gumbel"))
  )

  cells <- list()
  for (cfg in cfgs) {
    cat(sprintf("  cell %-22s ... ", cfg$name)); flush.console()
    t0   <- proc.time()
    cell <- run_mn_cell(cfg$n, cfg$betas, cfg$dgp_args, reps,
                        seed = 31415L + 7L * cfg$n + nchar(cfg$name))
    el   <- (proc.time() - t0)[3L]
    s    <- summarise_mn(cell)
    cat(sprintf("%5.1fs | df=%d  T-I(.05)=%.3f  mean=%.2f  KSp=%.3f\n",
                el, s$df, s$type1_05, s$mean_stat, s$ks_p))
    cells[[cfg$name]] <- list(meta = cfg, raw = cell, summary = s)
  }
  saveRDS(cells, file = out_path)
  cat("\nSaved to ", out_path, "\n", sep = "")
  invisible(cells)
}

if (sys.nframe() == 0L) {
  main_multinomial()
}
