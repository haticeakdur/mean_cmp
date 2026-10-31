# ============================================================
# FINAL REVISION SCRIPT 1: SIMULATION
# CMP_mu RQR/AQR calibration + dispersion score test
# ============================================================

required <- c("mpcmp", "nortest")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))

source("cmp_revision_helpers.R")

OUTDIR <- "revision_simulation_outputs"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Reproducibility
SEED_RQR   <- 123L
SEED_SCORE <- 459L
B_RQR      <- 2000L
B_SCORE    <- 2000L
ALPHA      <- 0.05

# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------
residual_metrics <- function(r) {
  r <- r[is.finite(r)]
  if (length(r) < 8L) return(rep(NA_real_, 5L))

  m <- mean(r)
  centered <- r - m
  m2 <- mean(centered^2)
  m3 <- mean(centered^3)
  m4 <- mean(centered^4)

  c(
    mean = m,
    variance = stats::var(r),
    skewness = if (m2 > 0) m3 / m2^(3/2) else NA_real_,
    excess_kurtosis = if (m2 > 0) m4 / m2^2 - 3 else NA_real_,
    AD = unname(nortest::ad.test(r)$statistic)
  )
}

# ------------------------------------------------------------
# Simulation 1: RQR and AQR under correctly specified CMP_mu
# ------------------------------------------------------------
run_rqr_aqr_scenario <- function(n, nu, B, scenario_id) {
  set.seed(SEED_RQR + 1000L * scenario_id)

  beta <- c(0.5, 0.7)
  store_rqr <- matrix(NA_real_, nrow = B, ncol = 5L)
  store_aqr <- matrix(NA_real_, nrow = B, ncol = 5L)
  colnames(store_rqr) <- colnames(store_aqr) <-
    c("mean", "variance", "skewness", "excess_kurtosis", "AD")

  n_fail <- 0L
  first_audit_error <- NA_real_
  first_generation_K <- NA_integer_

  for (b in seq_len(B)) {
    x <- stats::rnorm(n)
    mu <- exp(beta[1] + beta[2] * x)

    gen <- try(
      cmp_generate_mu(mu, nu, audit = (b == 1L)),
      silent = TRUE
    )
    if (inherits(gen, "try-error")) {
      n_fail <- n_fail + 1L
      next
    }
    if (b == 1L) {
      first_audit_error <- gen$audit_error
      first_generation_K <- gen$summax
    }

    dat <- data.frame(y = gen$y, x = x)
    fit <- try(
      mpcmp::glm.cmp(
        y ~ x,
        data = dat,
        lambdalb = 1e-10,
        lambdaub = 1000,
        tol = CMP_TOL,
        maxlambdaiter = CMP_MAX_LAMBDA_ITER
      ),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      n_fail <- n_fail + 1L
      next
    }

    # One randomized residual realization per fit; AQR uses the SAME RQR.
    rQ <- try(cmp_rqr(fit), silent = TRUE)
    if (inherits(rQ, "try-error")) {
      n_fail <- n_fail + 1L
      next
    }
    rAQ <- try(cmp_aqr(fit, rqr = rQ), silent = TRUE)
    if (inherits(rAQ, "try-error")) {
      n_fail <- n_fail + 1L
      next
    }

    store_rqr[b, ] <- residual_metrics(rQ)
    store_aqr[b, ] <- residual_metrics(rAQ)
  }

  summarize_one <- function(M, label) {
    ok <- stats::complete.cases(M)
    data.frame(
      Residual = label,
      n = n,
      nu = nu,
      Mean = mean(M[ok, "mean"]),
      Variance = mean(M[ok, "variance"]),
      Skewness = mean(M[ok, "skewness"]),
      Excess_Kurtosis = mean(M[ok, "excess_kurtosis"]),
      Mean_AD = mean(M[ok, "AD"]),
      SD_AD = stats::sd(M[ok, "AD"]),
      Successful = sum(ok),
      Failed = B - sum(ok),
      Generation_summax_first_rep = first_generation_K,
      Lambda_mean_constraint_max_error_first_rep = first_audit_error,
      stringsAsFactors = FALSE
    )
  }

  rbind(
    summarize_one(store_aqr, "AQR"),
    summarize_one(store_rqr, "RQR")
  )
}

rqr_grid <- expand.grid(
  n = c(50L, 200L),
  nu = c(0.6, 1.6),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

rqr_results <- do.call(
  rbind,
  lapply(seq_len(nrow(rqr_grid)), function(j) {
    cat("RQR/AQR scenario", j, "of", nrow(rqr_grid),
        ": n=", rqr_grid$n[j], ", nu=", rqr_grid$nu[j], "\n", sep = "")
    run_rqr_aqr_scenario(
      n = rqr_grid$n[j],
      nu = rqr_grid$nu[j],
      B = B_RQR,
      scenario_id = j
    )
  })
)

write.csv(
  rqr_results,
  file.path(OUTDIR, "Table1_RQR_AQR_revision.csv"),
  row.names = FALSE
)
print(rqr_results)

# ------------------------------------------------------------
# Simulation 2: dispersion score test
# H0: nu=1 vs H1: nu != 1
# ------------------------------------------------------------
run_score_scenario <- function(n, nu, B, scenario_id) {
  set.seed(SEED_SCORE + 1000L * scenario_id)

  beta0 <- 1
  beta1 <- 0.5
  reject <- rep(NA, B)
  first_audit_error <- NA_real_
  first_generation_K <- NA_integer_

  for (b in seq_len(B)) {
    x <- stats::runif(n, 0, 1)
    mu <- exp(beta0 + beta1 * x)

    gen <- try(
      cmp_generate_mu(mu, nu, audit = (b == 1L)),
      silent = TRUE
    )
    if (inherits(gen, "try-error")) next
    if (b == 1L) {
      first_audit_error <- gen$audit_error
      first_generation_K <- gen$summax
    }

    dat <- data.frame(y = gen$y, x = x)

    # Restricted MLE under H0: ordinary Poisson regression.
    fit0 <- try(stats::glm(y ~ x, data = dat, family = stats::poisson()), silent = TRUE)
    if (inherits(fit0, "try-error")) next

    mu0 <- stats::fitted(fit0)
    test <- try(cmp_score_test_null(y = dat$y, mu = mu0), silent = TRUE)
    if (inherits(test, "try-error")) next

    reject[b] <- test$p.value < ALPHA
  }

  data.frame(
    n = n,
    nu = nu,
    Dispersion = if (nu < 1) "Overdispersion" else if (nu > 1) "Underdispersion" else "Poisson (H0)",
    Rejection_probability = mean(reject, na.rm = TRUE),
    Successful = sum(!is.na(reject)),
    Failed = sum(is.na(reject)),
    Generation_summax_first_rep = first_generation_K,
    Lambda_mean_constraint_max_error_first_rep = first_audit_error,
    stringsAsFactors = FALSE
  )
}

score_grid <- expand.grid(
  n = c(100L, 200L),
  nu = c(1.6, 1.4, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

score_results <- do.call(
  rbind,
  lapply(seq_len(nrow(score_grid)), function(j) {
    cat("Score-test scenario", j, "of", nrow(score_grid),
        ": n=", score_grid$n[j], ", nu=", score_grid$nu[j], "\n", sep = "")
    run_score_scenario(
      n = score_grid$n[j],
      nu = score_grid$nu[j],
      B = B_SCORE,
      scenario_id = j
    )
  })
)

write.csv(
  score_results,
  file.path(OUTDIR, "Table2_score_test_revision.csv"),
  row.names = FALSE
)
print(score_results)

# Record exact software environment for the revision.
capture.output(
  sessionInfo(),
  file = file.path(OUTDIR, "sessionInfo_simulation.txt")
)

cat("\nSimulation revision completed. Outputs written to:", OUTDIR, "\n")
