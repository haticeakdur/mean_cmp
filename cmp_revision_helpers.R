# ============================================================
# Common helper functions for the CMP_mu revision
# Residual and Score-Based Diagnostics for Mean-Parametrized CMP Regression
# ============================================================

CMP_TOL <- 1e-6
CMP_MAX_LAMBDA_ITER <- 1000L

cmp_choose_summax <- function(mu, nu, min_k = 100L, mult = 20) {
  if (any(!is.finite(mu)) || any(mu <= 0)) stop("mu must be finite and > 0.")
  if (length(nu) != 1L || !is.finite(nu) || nu <= 0) stop("nu must be a finite positive scalar.")
  as.integer(ceiling(max(c(min_k, mu + mult * sqrt(mu / nu)))))
}

cmp_fit_summax <- function(object) {
  K <- object$summax
  if (is.null(K) || length(K) != 1L || !is.finite(K) || K < 2) {
    stop("The fitted CMP object does not contain a valid scalar summax.")
  }
  as.integer(K)
}

cmp_variance_fit <- function(object) {
  K <- cmp_fit_summax(object)
  V <- mpcmp::comp_variances(object$lambda, object$nu, object$log_Z, K)
  
  if (any(!is.finite(V)) || any(V <= 0)) {
    stop("Non-finite or non-positive CMP variance encountered; inspect numerical accuracy/truncation.")
  }
  
  V
}

rqr_from_bounds <- function(FL, FU, w = NULL, eps = .Machine$double.eps) {
  
  if (length(FL) != length(FU)) stop("FL and FU must have the same length.")
  
  if (any(!is.finite(FL)) || any(!is.finite(FU))) {
    stop("Non-finite CDF bounds encountered.")
  }
  
  if (any(FL < -1e-12) || any(FU > 1 + 1e-12)) {
    stop("CDF bounds fall outside [0,1].")
  }
  
  if (any(FU <= FL)) {
    stop("Invalid randomized-quantile interval: FU must be strictly greater than FL.")
  }
  
  n <- length(FL)
  
  if (is.null(w)) {
    w <- stats::runif(n)
  }
  
  if (length(w) != n ||
      any(!is.finite(w)) ||
      any(w < 0) ||
      any(w > 1)) {
    stop("w must be a vector in [0,1] with the same length as the CDF bounds.")
  }
  
  U <- FL + (FU - FL) * w
  
  # Numerical protection only after constructing
  # the exact randomized interval.
  U <- pmin(pmax(U, eps), 1 - eps)
  
  stats::qnorm(U)
}

# ============================================================
# Randomized quantile residual
# ============================================================

cmp_rqr <- function(object, w = NULL) {
  
  y <- object$y
  K <- cmp_fit_summax(object)
  
  FL <- mpcmp::pcomp(
    y - 1,
    lambda = object$lambda,
    nu = object$nu,
    summax = K
  )
  
  FU <- mpcmp::pcomp(
    y,
    lambda = object$lambda,
    nu = object$nu,
    summax = K
  )
  
  rqr_from_bounds(
    FL = FL,
    FU = FU,
    w = w
  )
}

# ============================================================
# Adjusted randomized quantile residual
# Same random realization of the RQR can be supplied
# ============================================================

cmp_aqr <- function(object, rqr = NULL) {
  
  if (is.null(rqr)) {
    rqr <- cmp_rqr(object)
  }
  
  h <- stats::hatvalues(object)
  
  if (length(h) != length(rqr)) {
    stop("Leverage and residual vectors have different lengths.")
  }
  
  if (any(!is.finite(h)) ||
      any(h < 0) ||
      any(h >= 1)) {
    stop("Invalid leverage values encountered in AQR calculation.")
  }
  
  rqr / sqrt(1 - h)
}

# ============================================================
# Pearson residual
# ============================================================

cmp_pearson <- function(object) {
  
  (object$y - object$fitted_values) /
    sqrt(cmp_variance_fit(object))
}

# ============================================================
# Leverage-adjusted Pearson residual
# This is NOT the Williams likelihood residual
# ============================================================

cmp_lap <- function(object) {
  
  rP <- cmp_pearson(object)
  h <- stats::hatvalues(object)
  
  if (any(!is.finite(h)) ||
      any(h < 0) ||
      any(h >= 1)) {
    stop("Invalid leverage values encountered in leverage-adjusted Pearson residual.")
  }
  
  rP / sqrt(1 - h)
}

# ============================================================
# Anscombe-type residual
#
# A_nu(mu) = integral V(t,nu)^(-1/3) dt
#
# r_A =
# [A_nu(y) - A_nu(mu)] /
# V(mu,nu)^(1/6)
#
# The integral is evaluated separately at fitted nu
# using a trapezoidal rule.
# ============================================================

cmp_anscombe <- function(object,
                         n_grid = 1200L,
                         eps = 1e-8,
                         tol = CMP_TOL,
                         maxlambdaiter = CMP_MAX_LAMBDA_ITER) {
  
  y  <- object$y
  mu <- object$fitted_values
  nu <- object$nu
  
  if (length(nu) != 1L ||
      !is.finite(nu) ||
      nu <= 0) {
    stop("cmp_anscombe currently assumes a constant positive dispersion parameter.")
  }
  
  if (n_grid < 200L) {
    stop("Use at least 200 grid points for numerical integration.")
  }
  
  upper <- 1.05 * max(c(mu, y, 1), na.rm = TRUE)
  
  u <- seq(
    from = 0,
    to = 1,
    length.out = n_grid
  )
  
  # Quadratic spacing places more points near zero,
  # where V(mu,nu)^(-1/3) may vary rapidly.
  grid <- eps + (upper - eps) * u^2
  
  # Include fitted means and observed counts
  # explicitly in the integration grid.
  grid <- sort(
    unique(
      c(
        grid,
        mu[mu > eps],
        y[y > eps]
      )
    )
  )
  
  K <- max(
    cmp_fit_summax(object),
    cmp_choose_summax(grid, nu)
  )
  
  # For every grid mean, solve lambda(mu,nu)
  # from the mean constraint.
  lambda_grid <- mpcmp::comp_lambdas(
    mu = grid,
    nu = nu,
    lambdalb = 1e-10,
    lambdaub = 1000,
    maxlambdaiter = maxlambdaiter,
    tol = tol,
    summax = K
  )$lambda
  
  if (any(!is.finite(lambda_grid)) ||
      any(lambda_grid <= 0)) {
    stop("Lambda inversion failed on the Anscombe integration grid.")
  }
  
  logZ_grid <- mpcmp::logZ(
    log(lambda_grid),
    nu,
    summax = K
  )
  
  V_grid <- mpcmp::comp_variances(
    lambda_grid,
    nu,
    logZ_grid,
    K
  )
  
  if (any(!is.finite(V_grid)) ||
      any(V_grid <= 0)) {
    stop("Invalid variances on the Anscombe integration grid.")
  }
  
  integrand <- V_grid^(-1 / 3)
  
  dmu <- diff(grid)
  
  # Near zero, V(mu,nu) ~ mu.
  # Therefore:
  #
  # integral_0^eps t^(-1/3) dt
  # = (3/2) eps^(2/3)
  #
  A_eps <- 1.5 * eps^(2 / 3)
  
  A_grid <- A_eps + c(
    0,
    cumsum(
      0.5 *
        (
          integrand[-1] +
            integrand[-length(integrand)]
        ) *
        dmu
    )
  )
  
  A_fun <- stats::approxfun(
    x = c(0, grid),
    y = c(0, A_grid),
    method = "linear",
    rule = 2,
    ties = "ordered"
  )
  
  V_obs <- cmp_variance_fit(object)
  
  rA <- (
    A_fun(y) -
      A_fun(mu)
  ) / V_obs^(1 / 6)
  
  if (any(!is.finite(rA))) {
    stop("Non-finite Anscombe-type residual encountered.")
  }
  
  rA
}

# ============================================================
# Dispersion-score residual
#
# Individual dispersion-score contributions:
#
# S_nu,i =
# A_i (y_i - mu_i) / V_i
# - [log(y_i!) - B_i]
#
# Empirical standardization is only for graphical
# comparability. It does NOT imply N(0,1).
# ============================================================

cmp_dispersion_score <- function(object,
                                 standardize = TRUE) {
  
  y      <- object$y
  lambda <- object$lambda
  nu     <- object$nu
  logZ   <- object$log_Z
  
  K <- cmp_fit_summax(object)
  
  mu_vec <- mpcmp::comp_means(
    lambda,
    nu,
    logZ,
    K
  )
  
  B_vec <- mpcmp::comp_mean_logfactorialy(
    lambda,
    nu,
    logZ,
    K
  )
  
  Eylog_vec <- mpcmp::comp_mean_ylogfactorialy(
    lambda,
    nu,
    logZ,
    K
  )
  
  V_vec <- mpcmp::comp_variances(
    lambda,
    nu,
    logZ,
    K
  )
  
  if (any(!is.finite(V_vec)) ||
      any(V_vec <= 0)) {
    stop("Invalid variance encountered in dispersion-score calculation.")
  }
  
  A_vec <- Eylog_vec -
    mu_vec * B_vec
  
  S_nu <- A_vec *
    (y - mu_vec) /
    V_vec -
    (
      lgamma(y + 1) -
        B_vec
    )
  
  if (any(!is.finite(S_nu))) {
    stop("Non-finite dispersion-score contribution encountered.")
  }
  
  if (!standardize) {
    return(S_nu)
  }
  
  S_sd <- stats::sd(S_nu)
  
  if (!is.finite(S_sd) ||
      S_sd <= 0) {
    stop("Dispersion-score contributions have zero/invalid empirical SD.")
  }
  
  (
    S_nu -
      mean(S_nu)
  ) / S_sd
}

# ============================================================
# Dispersion score under H0: nu = 1
#
# Under H0:
# lambda = mu
# V = mu
# log Z = mu
#
# B and E[Y log(Y!)] are evaluated using
# the same truncated mpcmp moment calculations.
# ============================================================

cmp_score_null_contrib <- function(y,
                                   mu,
                                   summax = NULL) {
  
  if (length(y) != length(mu)) {
    stop("y and mu must have the same length.")
  }
  
  if (any(!is.finite(mu)) ||
      any(mu <= 0)) {
    stop("Null fitted means must be finite and > 0.")
  }
  
  if (is.null(summax)) {
    summax <- cmp_choose_summax(
      mu = mu,
      nu = 1
    )
  }
  
  K <- as.integer(summax)
  
  # Poisson special case
  lambda0 <- mu
  
  # Z(lambda,1)=exp(lambda)
  # therefore log Z=lambda=mu
  logZ0 <- mu
  
  B_vec <- mpcmp::comp_mean_logfactorialy(
    lambda0,
    1,
    logZ0,
    K
  )
  
  Eylog_vec <- mpcmp::comp_mean_ylogfactorialy(
    lambda0,
    1,
    logZ0,
    K
  )
  
  A_vec <- Eylog_vec -
    mu * B_vec
  
  u_nu <- A_vec *
    (y - mu) /
    mu -
    (
      lgamma(y + 1) -
        B_vec
    )
  
  if (any(!is.finite(u_nu))) {
    stop("Non-finite null dispersion-score contribution encountered.")
  }
  
  u_nu
}

# ============================================================
# Dispersion score test
#
# SC =
# (sum u_i)^2 / sum u_i^2
#
# OPG estimator of information
# ============================================================

cmp_score_test_null <- function(y,
                                mu,
                                summax = NULL) {
  
  u <- cmp_score_null_contrib(
    y = y,
    mu = mu,
    summax = summax
  )
  
  U <- sum(u)
  
  Iopg <- sum(u^2)
  
  if (!is.finite(Iopg) ||
      Iopg <= 0) {
    stop("Invalid OPG information in score test.")
  }
  
  SC <- U^2 / Iopg
  
  list(
    statistic = unname(SC),
    df = 1L,
    p.value = stats::pchisq(
      SC,
      df = 1,
      lower.tail = FALSE
    ),
    score_sum = U,
    opg_information = Iopg,
    contributions = u,
    summax = if (is.null(summax)) {
      cmp_choose_summax(mu, 1)
    } else {
      as.integer(summax)
    }
  )
}

# ============================================================
# CMP_mu data generation
#
# mpcmp::rcomp(mu=..., nu=...) is used.
#
# The optional audit checks the transformation:
#
# mu -> lambda(mu,nu) -> E(Y | lambda,nu)
#
# and verifies that the recovered mean is numerically
# equal to the target mu.
# ============================================================

cmp_generate_mu <- function(mu,
                            nu,
                            tol = CMP_TOL,
                            maxlambdaiter = CMP_MAX_LAMBDA_ITER,
                            summax = NULL,
                            audit = FALSE) {
  
  if (is.null(summax)) {
    summax <- cmp_choose_summax(
      mu = mu,
      nu = nu
    )
  }
  
  K <- as.integer(summax)
  
  audit_error <- NA_real_
  
  if (audit) {
    
    solution <- mpcmp::comp_lambdas(
      mu = mu,
      nu = nu,
      lambdalb = 1e-10,
      lambdaub = 1000,
      maxlambdaiter = maxlambdaiter,
      tol = tol,
      summax = K
    )
    
    lambda <- solution$lambda
    
    logZ <- mpcmp::logZ(
      log(lambda),
      nu,
      summax = K
    )
    
    mu_check <- mpcmp::comp_means(
      lambda,
      nu,
      logZ,
      K
    )
    
    audit_error <- max(
      abs(mu_check - mu)
    )
    
    if (!is.finite(audit_error) ||
        audit_error > 5 * tol) {
      stop(
        sprintf(
          "Mean-to-lambda inversion audit failed: max error = %.3e",
          audit_error
        )
      )
    }
  }
  
  y <- mpcmp::rcomp(
    n = length(mu),
    mu = mu,
    nu = nu,
    lambdalb = 1e-10,
    lambdaub = 1000,
    maxlambdaiter = maxlambdaiter,
    tol = tol,
    summax = K
  )
  
  if (any(!is.finite(y))) {
    stop("CMP_mu generation returned non-finite values.")
  }
  
  list(
    y = as.integer(y),
    summax = K,
    audit_error = audit_error
  )
}