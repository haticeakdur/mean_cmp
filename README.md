# Residual and Score-Based Diagnostics for Mean-Parametrized Conway–Maxwell–Poisson Regression

This repository contains the reproducibility code for the manuscript **“Residual and Score-Based Diagnostics for Mean-Parametrized Conway–Maxwell–Poisson Regression.”**

The code is organized so that the simulation study and empirical application can be reproduced directly from the repository root.

## Files

- `cmp_revision_helpers.R`  
  Common CMP_mu numerical and diagnostic functions used by both analyses. This file contains:
  - the numerical tolerance and maximum number of lambda-iterations;
  - the truncation rule for CMP sums;
  - mean-to-rate inversion through `mpcmp::comp_lambdas()`;
  - CMP variance calculation through `mpcmp::comp_variances()`;
  - CMP CDF evaluation through `mpcmp::pcomp()`;
  - randomized quantile residual (RQR) and adjusted RQR (AQR);
  - Pearson and leverage-adjusted Pearson (LAP) residuals;
  - the numerical Anscombe-type residual;
  - individual dispersion-score contributions and the dispersion score test;
  - mean-parametrized CMP data generation and an optional numerical audit of the mean constraint.

- `cmp_simulation_revision.R`  
  Exact code used for the Monte Carlo studies reported in the revised manuscript.

- `cmp_application_revision.R`  
  Exact code used for the takeover-bids application, model comparison, RQR calibration, CMP_mu-specific diagnostics, exact single-case deletion analysis, and the dispersion score test.

- `revision_simulation_outputs/`  
  Created by the simulation script. It contains the revised simulation tables and `sessionInfo_simulation.txt`.

- `revision_application_outputs/`  
  Contains the application tables, diagnostic data, figures, influence summaries, score-test output, and `sessionInfo_application.txt`.

## Reproducibility settings

### Random seeds

Simulation study:

- RQR/AQR study base seed: `123`
- dispersion score-test study base seed: `459`
- within each simulation design point, the script uses `base seed + 1000 * scenario_id`

Empirical application:

- RQR seed: `20261031`

### Numerical settings for CMP_mu

The numerical settings used throughout the revised analyses are defined explicitly in `cmp_revision_helpers.R`:

- lower lambda search bound: `1e-10`
- upper lambda search bound: `1000`
- convergence tolerance: `1e-6`
- maximum lambda iterations: `1000`

The same bounds, tolerance, and maximum iteration count are specified explicitly in the simulation and empirical `glm.cmp()` fits.

### Truncation rule for CMP sums

When an explicit truncation level is needed, the helper function uses

```r
K = ceiling(max(100, mu + 20 * sqrt(mu / nu)))
```

with the maximum taken over the relevant mean vector. For fitted CMP_mu objects, the truncation level stored in the fitted model is used for CDF and fitted-moment calculations. For the Anscombe integration grid, the truncation level is chosen to be at least as large as the fitted-model value.

### Mean-to-rate inversion

For mean-parametrized CMP calculations, the conventional CMP rate parameter `lambda(mu, nu)` is obtained numerically by solving

```text
E(Y | lambda, nu) = mu
```

using `mpcmp::comp_lambdas()` with the numerical settings listed above. In the simulation data-generation routine, `mpcmp::rcomp(mu = ..., nu = ...)` is used; for the first replication of each design point, an additional audit recovers `lambda` using `comp_lambdas()` and checks the resulting CMP mean against the target `mu`.

### Variance function

The CMP_mu variance is evaluated from truncated CMP moments using

```r
mpcmp::comp_variances()
```

with the fitted `lambda`, `nu`, normalizing constant, and truncation level.

### CDF and randomized quantile residuals

For the CMP_mu model, the two CDF bounds used in the randomized probability integral transform are evaluated as

```r
F(y_i - 1; lambda_i, nu)
F(y_i;     lambda_i, nu)
```

using `mpcmp::pcomp()`. The same randomized realization is used for the RQR and its leverage-adjusted version AQR.

For Poisson, NB2, generalized Poisson, and Bell models in the empirical application, Dunn–Smyth residuals are obtained directly from `glmmTMB` using

```r
residuals(fit, type = "dunn-smyth")
```

No custom Bell or generalized-Poisson RQR implementation is used in the revised analysis.

## Software environment

The empirical analysis was run under:

- R 4.5.1
- `mpcmp` 0.3.8
- `glmmTMB` 1.1.15
- `ggplot2` 4.0.1
- `nortest` 1.0-4
- `patchwork` 1.3.2

The complete software environment is recorded in `revision_application_outputs/sessionInfo_application.txt`. The simulation script similarly writes `revision_simulation_outputs/sessionInfo_simulation.txt` when run.

## Running the analyses

From the repository root, run:

```r
source("cmp_simulation_revision.R")
source("cmp_application_revision.R")
```

Both scripts source `cmp_revision_helpers.R` automatically and create their corresponding output directories.

## Data

The takeover-bids dataset is loaded directly from the `mpcmp` package:

```r
data("takeoverbids", package = "mpcmp")
```

No external data download is required.

## Repository link

https://github.com/haticeakdur/mean_cmp
