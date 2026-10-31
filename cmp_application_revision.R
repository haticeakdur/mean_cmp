# ============================================================
# FINAL REVISION SCRIPT 2: EMPIRICAL APPLICATION
# Takeover bids data
# ============================================================

required <- c(
  "mpcmp", "glmmTMB", "ggplot2", "nortest", "patchwork"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))

source("cmp_revision_helpers.R")

OUTDIR <- "revision_application_outputs"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Data and common model specification
# ------------------------------------------------------------
data("takeoverbids", package = "mpcmp")

vars <- c(
  "numbids", "leglrest", "rearest", "finrest", "whtknght",
  "bidprem", "insthold", "size", "sizesq", "regulatn"
)

dat <- takeoverbids[stats::complete.cases(takeoverbids[, vars]), vars]
rownames(dat) <- NULL

form_bids <- numbids ~ leglrest + rearest + finrest + whtknght +
  bidprem + insthold + size + sizesq + regulatn

# All models are fitted to EXACTLY the same rows.
fit_pois <- glmmTMB::glmmTMB(form_bids, data = dat, family = stats::poisson())
fit_nb2  <- glmmTMB::glmmTMB(form_bids, data = dat, family = glmmTMB::nbinom2())
fit_gp   <- glmmTMB::glmmTMB(form_bids, data = dat, family = glmmTMB::genpois())
fit_bell <- glmmTMB::glmmTMB(form_bids, data = dat, family = glmmTMB::bell())
fit_cmp  <- mpcmp::glm.cmp(
  form_bids,
  data = dat,
  lambdalb = 1e-10,
  lambdaub = 1000,
  tol = CMP_TOL,
  maxlambdaiter = CMP_MAX_LAMBDA_ITER
)

extract_fit_stats <- function(fit, model_name) {
  ll <- stats::logLik(fit)
  data.frame(
    Model = model_name,
    df = attr(ll, "df"),
    LogLik = as.numeric(ll),
    AIC = stats::AIC(fit),
    stringsAsFactors = FALSE
  )
}

model_table <- rbind(
  extract_fit_stats(fit_pois, "Poisson"),
  extract_fit_stats(fit_nb2,  "NB2"),
  extract_fit_stats(fit_gp,   "Generalized Poisson"),
  extract_fit_stats(fit_bell, "Bell"),
  extract_fit_stats(fit_cmp,  "CMP_mu")
)
model_table$Delta_AIC <- model_table$AIC - min(model_table$AIC)
model_table <- model_table[order(model_table$AIC), ]
rownames(model_table) <- NULL

write.csv(model_table, file.path(OUTDIR, "application_model_comparison.csv"), row.names = FALSE)
print(model_table)

nu_hat <- fit_cmp$nu
nu_interpretation <- if (nu_hat > 1) "underdispersion" else if (nu_hat < 1) "overdispersion" else "equidispersion"
cat(sprintf("\nCMP_mu nu_hat = %.6f -> %s.\n", nu_hat, nu_interpretation))

# ------------------------------------------------------------
# 2. Randomized quantile (Dunn-Smyth) residuals
# ------------------------------------------------------------
# Current glmmTMB implements type = "dunn-smyth" directly for
# Poisson, NB2, generalized Poisson, and Bell families.
# The same seed is reset before each calculation so that the
# randomization is reproducible and based on the same RNG stream.
RQR_SEED <- 20261031L
y <- dat$numbids

set.seed(RQR_SEED)
rqr_pois <- stats::residuals(fit_pois, type = "dunn-smyth")

set.seed(RQR_SEED)
rqr_nb2 <- stats::residuals(fit_nb2, type = "dunn-smyth")

set.seed(RQR_SEED)
rqr_gp <- stats::residuals(fit_gp, type = "dunn-smyth")

set.seed(RQR_SEED)
rqr_bell <- stats::residuals(fit_bell, type = "dunn-smyth")

# CMP_mu RQR: evaluated from the fitted lambda_i, nu, and summax
# stored by glm.cmp, using the same underlying Uniform(0,1) stream.
set.seed(RQR_SEED)
w_cmp <- stats::runif(nrow(dat))
rqr_cmp <- cmp_rqr(fit_cmp, w = w_cmp)

rqr_df <- data.frame(
  Model = rep(c("Poisson", "NB2", "Generalized Poisson", "Bell", "CMP_mu"), each = nrow(dat)),
  RQR = c(rqr_pois, rqr_nb2, rqr_gp, rqr_bell, rqr_cmp),
  stringsAsFactors = FALSE
)

calibration_one <- function(r, label) {
  r <- r[is.finite(r)]
  ad <- nortest::ad.test(r)
  data.frame(
    Model = label,
    Mean = mean(r),
    Variance = stats::var(r),
    AD_statistic = unname(ad$statistic),
    AD_p_value_descriptive_only = ad$p.value,
    stringsAsFactors = FALSE
  )
}

rqr_calibration <- rbind(
  calibration_one(rqr_pois, "Poisson"),
  calibration_one(rqr_nb2, "NB2"),
  calibration_one(rqr_gp, "Generalized Poisson"),
  calibration_one(rqr_bell, "Bell"),
  calibration_one(rqr_cmp, "CMP_mu")
)
write.csv(rqr_calibration, file.path(OUTDIR, "application_RQR_calibration.csv"), row.names = FALSE)
print(rqr_calibration)

p_qq <- ggplot2::ggplot(rqr_df, ggplot2::aes(sample = RQR)) +
  ggplot2::stat_qq(alpha = 0.6) +
  ggplot2::stat_qq_line(linetype = "dashed") +
  ggplot2::facet_wrap(~ Model) +
  ggplot2::labs(
    title = "Randomized Quantile Residuals Across Count Models",
    x = "Theoretical quantiles",
    y = "Sample quantiles"
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(file.path(OUTDIR, "Figure1_RQR_QQ_revision.pdf"), p_qq, width = 10, height = 6)

# ------------------------------------------------------------
# 3. CMP_mu-specific residual diagnostics
# ------------------------------------------------------------
rP   <- cmp_pearson(fit_cmp)
rLAP <- cmp_lap(fit_cmp)
rA   <- cmp_anscombe(fit_cmp)
ds_raw <- cmp_dispersion_score(fit_cmp, standardize = FALSE)
ds_std <- cmp_dispersion_score(fit_cmp, standardize = TRUE)
h <- stats::hatvalues(fit_cmp)
mu_cmp <- fit_cmp$fitted_values

cmp_diag <- data.frame(
  obs = seq_len(nrow(dat)),
  y = y,
  mu_hat = mu_cmp,
  leverage = h,
  Pearson = rP,
  LAP = rLAP,
  Anscombe = rA,
  DispersionScore_raw = ds_raw,
  DispersionScore_std = ds_std
)
write.csv(cmp_diag, file.path(OUTDIR, "application_CMP_diagnostics.csv"), row.names = FALSE)

# Targeted diagnostic panels
p_lap <- ggplot2::ggplot(cmp_diag, ggplot2::aes(leverage, LAP)) +
  ggplot2::geom_point(alpha = 0.65) +
  ggplot2::geom_smooth(method = "loess", se = FALSE) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::labs(title = "Leverage-adjusted Pearson residual", x = "Leverage", y = "LAP residual") +
  ggplot2::theme_minimal(base_size = 12)

p_ans <- ggplot2::ggplot(cmp_diag, ggplot2::aes(mu_hat, Anscombe)) +
  ggplot2::geom_point(alpha = 0.65) +
  ggplot2::geom_smooth(method = "loess", se = FALSE) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::labs(title = "Anscombe-type residual", x = expression(hat(mu)), y = "Anscombe residual") +
  ggplot2::theme_minimal(base_size = 12)

p_ds <- ggplot2::ggplot(cmp_diag, ggplot2::aes(mu_hat, DispersionScore_std)) +
  ggplot2::geom_point(alpha = 0.65) +
  ggplot2::geom_smooth(method = "loess", se = FALSE) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::labs(
    title = "Dispersion-score contribution",
    x = expression(hat(mu)),
    y = "Standardized score contribution"
  ) +
  ggplot2::theme_minimal(base_size = 12)

p_diag <- p_lap + p_ans + p_ds + patchwork::plot_layout(ncol = 3)
ggplot2::ggsave(file.path(OUTDIR, "Figure2_CMP_targeted_diagnostics_revision.pdf"), p_diag, width = 14, height = 4.8)

# ------------------------------------------------------------
# 4. Influence analysis: exact single-case deletion
#    NO AIC comparison after deletion.
# ------------------------------------------------------------
loo_cmp_influence <- function(fit_full, data, formula) {
  n <- nrow(data)
  beta_full <- stats::coef(fit_full)
  Vbeta <- stats::vcov(fit_full)
  if (any(!is.finite(Vbeta))) stop("Invalid covariance matrix for full CMP fit.")
  Vbeta_inv <- solve(Vbeta)
  se_beta <- sqrt(diag(Vbeta))
  p <- length(beta_full)
  nu_full <- fit_full$nu
  mu_full <- fit_full$fitted_values

  out <- vector("list", n)

  for (i in seq_len(n)) {
    cat("Leave-one-out CMP fit", i, "of", n, "\r")
    fit_i <- try(
      mpcmp::glm.cmp(
        formula,
        data = data[-i, , drop = FALSE],
        lambdalb = 1e-10,
        lambdaub = 1000,
        tol = CMP_TOL,
        maxlambdaiter = CMP_MAX_LAMBDA_ITER
      ),
      silent = TRUE
    )

    if (inherits(fit_i, "try-error")) {
      out[[i]] <- data.frame(
        obs = i, Cook_type_beta = NA_real_, Max_abs_std_beta_change = NA_real_,
        Delta_nu = NA_real_, Relative_nu_change_percent = NA_real_,
        RMS_fitted_change = NA_real_, Max_abs_fitted_change = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    beta_i <- stats::coef(fit_i)[names(beta_full)]
    delta_beta <- beta_i - beta_full
    pred_i <- as.numeric(stats::predict(fit_i, newdata = data, type = "response"))

    cook_type_beta <- as.numeric(
      t(delta_beta) %*% Vbeta_inv %*% delta_beta / p
    )

    out[[i]] <- data.frame(
      obs = i,
      Cook_type_beta = cook_type_beta,
      Max_abs_std_beta_change = max(abs(delta_beta / se_beta)),
      Delta_nu = fit_i$nu - nu_full,
      Relative_nu_change_percent = 100 * (fit_i$nu - nu_full) / nu_full,
      RMS_fitted_change = sqrt(mean((pred_i - mu_full)^2)),
      Max_abs_fitted_change = max(abs(pred_i - mu_full)),
      stringsAsFactors = FALSE
    )
  }
  cat("\n")
  do.call(rbind, out)
}

influence_tbl <- loo_cmp_influence(fit_cmp, dat, form_bids)
influence_tbl$y <- y[influence_tbl$obs]
influence_tbl$leverage <- h[influence_tbl$obs]
influence_tbl$LAP <- rLAP[influence_tbl$obs]
influence_tbl$DispersionScore_std <- ds_std[influence_tbl$obs]

write.csv(influence_tbl, file.path(OUTDIR, "application_LOO_influence_all.csv"), row.names = FALSE)

top10_beta <- influence_tbl[order(influence_tbl$Cook_type_beta, decreasing = TRUE), ][1:10, ]
top10_nu <- influence_tbl[order(abs(influence_tbl$Delta_nu), decreasing = TRUE), ][1:10, ]

write.csv(top10_beta, file.path(OUTDIR, "application_LOO_top10_beta_influence.csv"), row.names = FALSE)
write.csv(top10_nu, file.path(OUTDIR, "application_LOO_top10_nu_influence.csv"), row.names = FALSE)

cat("\nTop 10 observations by Cook-type coefficient distance:\n")
print(top10_beta)
cat("\nTop 10 observations by |Delta nu|:\n")
print(top10_nu)

p_inf <- ggplot2::ggplot(influence_tbl, ggplot2::aes(obs, Cook_type_beta)) +
  ggplot2::geom_point(alpha = 0.7) +
  ggplot2::labs(
    title = "Exact single-case deletion influence for CMP_mu regression",
    x = "Observation index",
    y = "Cook-type coefficient distance"
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(file.path(OUTDIR, "CMP_LOO_Cook_type_influence.pdf"), p_inf, width = 8, height = 5)

# ------------------------------------------------------------
# 5. Dispersion score test for H0: nu=1
# ------------------------------------------------------------
fit_pois_null <- stats::glm(form_bids, data = dat, family = stats::poisson())
mu0 <- stats::fitted(fit_pois_null)
score_test <- cmp_score_test_null(y = dat$numbids, mu = mu0)

score_result <- data.frame(
  statistic = score_test$statistic,
  df = score_test$df,
  p_value = score_test$p.value,
  summax = score_test$summax,
  stringsAsFactors = FALSE
)
write.csv(score_result, file.path(OUTDIR, "application_score_test.csv"), row.names = FALSE)
print(score_result)

# ------------------------------------------------------------
# 6. Reproducibility record
# ------------------------------------------------------------
capture.output(
  sessionInfo(),
  file = file.path(OUTDIR, "sessionInfo_application.txt")
)

cat("\nApplication revision completed. Outputs written to:", OUTDIR, "\n")
