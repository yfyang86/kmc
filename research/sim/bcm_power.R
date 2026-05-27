## sim/bcm_power.R
##
## Power study for the BCM Owen-EL test (calibration in bcm_sim.R).
## Same machinery; here we generate data from the truth beta_true but
## test H0: beta = beta_true + (0, h, 0, ...) for a sequence of h.

source("sim/bcm_sim.R", chdir = FALSE)

run_power_cell <- function(n, beta_true, dgp, h_grid, reps, seed) {
  set.seed(seed)
  power <- numeric(length(h_grid))
  crit  <- qchisq(0.95, df = length(beta_true) - 1L)
  for (k in seq_along(h_grid)) {
    h     <- h_grid[k]
    rej   <- 0L
    for (r in seq_len(reps)) {
      d <- dgp(n, beta_true)
      # Test H0: beta = beta_true + h * e_2  (perturb the SECOND component)
      beta0 <- beta_true
      beta0[2L] <- beta0[2L] + h
      out <- bcm_lr_owen(d$X, d$delta, beta0)
      if (!is.na(out$stat) && out$stat > crit) rej <- rej + 1L
    }
    power[k] <- rej / reps
  }
  data.frame(h = h_grid, power = power)
}

main_power <- function(reps     = 2000L,
                       out_path = "sim/results/power_results.rds") {
  cat(sprintf("\nBCM Owen-EL power study: reps=%d per (cell, h)\n", reps))
  cat(strrep("=", 70), "\n", sep = "")

  beta_p2 <- c(1.0,  0.5)
  beta_p3 <- c(1.0,  0.5, -0.3)
  h_grid  <- c(-0.50, -0.30, -0.15, -0.05, 0.0,
                0.05, 0.15, 0.30, 0.50, 0.75, 1.00)

  cfgs <- list(
    list(name = "logistic_p2_n500",  n = 500L,  beta = beta_p2, dgp = gen_logistic),
    list(name = "logistic_p2_n1000", n = 1000L, beta = beta_p2, dgp = gen_logistic),
    list(name = "logistic_p2_n2000", n = 2000L, beta = beta_p2, dgp = gen_logistic),
    list(name = "logistic_p3_n1000", n = 1000L, beta = beta_p3, dgp = gen_logistic)
  )

  cells <- list()
  for (cfg in cfgs) {
    cat(sprintf("  cell %-22s ... ", cfg$name)); flush.console()
    t0 <- proc.time()
    pw <- run_power_cell(cfg$n, cfg$beta, cfg$dgp, h_grid, reps,
                         seed = 3026L + 13L * cfg$n + nchar(cfg$name))
    el <- (proc.time() - t0)[3L]
    cat(sprintf("%6.1fs\n", el))
    print(pw)
    cells[[cfg$name]] <- list(meta = cfg, power = pw)
  }

  saveRDS(cells, file = out_path)
  cat("\nSaved to ", out_path, "\n", sep = "")
  invisible(cells)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  reps <- if (length(args) >= 1L) as.integer(args[1L]) else 2000L
  outp <- if (length(args) >= 2L) args[2L] else "sim/results/power_results.rds"
  main_power(reps = reps, out_path = outp)
}
