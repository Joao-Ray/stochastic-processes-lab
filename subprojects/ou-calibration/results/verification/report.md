# OU Numerical Verification: Executed Report

Derivations and assumptions: [THEORY.md](../../THEORY.md). Full results: [results.json](results.json).

Seed 20261009; 50,000 paths; eight independently seeded likelihood records.

## Itô variance versus adaptive quadrature

| theta | h | Formula | Quadrature | Absolute error |
| --- | --- | --- | --- | --- |
| 0.05 | 0.01 | 0.006396801 | 0.006396801 | 1.93e-16 |
| 0.05 | 0.2 | 0.126728491 | 0.126728491 | 3.05e-16 |
| 0.05 | 2.0 | 1.160123180 | 1.160123180 | 4.44e-16 |
| 0.7 | 0.01 | 0.006355408 | 0.006355408 | 6.07e-18 |
| 0.7 | 0.2 | 0.111641718 | 0.111641718 | 0.00e+00 |
| 0.7 | 2.0 | 0.429343971 | 0.429343971 | 5.55e-17 |
| 3.0 | 0.01 | 0.006211783 | 0.006211783 | 1.73e-18 |
| 3.0 | 0.2 | 0.074539284 | 0.074539284 | 0.00e+00 |
| 3.0 | 2.0 | 0.106666011 | 0.106666011 | 1.39e-17 |

## Conditional moments and known-parameter prediction intervals

| h | Mean theory / empirical | Variance theory / empirical | Mean error / SE | Variance error / SE | Coverage [Wilson 95%] |
| --- | --- | --- | --- | --- | --- |
| 0.1 | -0.830985 / -0.831211 | 0.059722 / 0.059392 | -0.207 | -0.872 | 0.9509 [0.9489, 0.9527] |
| 0.5 | -0.261720 / -0.264258 | 0.230132 / 0.231310 | -1.183 | 0.809 | 0.9497 [0.9478, 0.9516] |
| 1.0 | 0.258537 / 0.249572 | 0.344413 / 0.346134 | -3.416 | 0.790 | 0.9495 [0.9475, 0.9513] |
| 3.0 | 1.193859 / 1.193678 | 0.450288 / 0.451659 | -0.060 | 0.481 | 0.9497 [0.9477, 0.9516] |

Moments condition on the fixed initial value x0. Horizons share paths, so errors across rows are correlated. These are known-parameter intervals, not fitted plug-in intervals.

Pointwise 95% mean intervals miss their theoretical target at 1 of 4 horizons: h=1.0 (-3.42 SE). These outcomes are retained; shared horizons are not independent tests.

## Conditional MSE: optimal prediction versus persistence

| h | Oracle theory / empirical | Persistence theory / empirical |
| --- | --- | --- |
| 0.1 | 0.059722 / 0.059391 | 0.088288 / 0.087881 |
| 0.5 | 0.230132 / 0.231311 | 0.775189 / 0.772621 |
| 1.0 | 0.344413 / 0.346208 | 1.928328 / 1.907558 |
| 3.0 | 0.450288 / 0.451650 | 5.263305 / 5.263872 |

![MSE verification and moment sampling errors](forecast_mse_verification.png)

## Conditional MLE: closed form versus numerical minimization

Each numerical reference uses two fixed starts in log(theta), mu, log(sigma), with wide bounds recorded in JSON. No start uses the closed-form estimate. Optimization statuses for both starts are retained; agreement is judged by objective and parameter errors, not the status flag alone.

| Seed | Absolute NLL gap | Max scaled parameter error | Normal equations error | Start convergence flags |
| --- | --- | --- | --- | --- |
| 20261109 | 0.00e+00 | 1.66e-08 | 1.91e-16 | True, True |
| 20261110 | 2.27e-13 | 2.95e-08 | 3.23e-16 | True, True |
| 20261111 | 4.55e-13 | 8.17e-08 | 2.92e-16 | True, True |
| 20261112 | 5.68e-14 | 6.83e-09 | 2.40e-16 | True, True |
| 20261113 | 1.71e-13 | 5.35e-08 | 2.46e-16 | True, True |
| 20261114 | 1.71e-13 | 4.46e-08 | 2.70e-16 | True, True |
| 20261115 | 5.68e-14 | 1.26e-08 | 2.10e-16 | True, True |
| 20261116 | 5.68e-13 | 4.79e-08 | 1.15e-16 | True, True |

## Composition and paired mean shift

- Maximum 50-step recursive-moment discrepancy: **8.88e-16**.
- Maximum coupled-path mean-shift identity error: **4.66e-15**.
- Numerical ODE versus closed-form displacement error: **2.57e-11**.
- Old known-parameter interval coverage after shift: analytic **0.4236**, simulated **0.4258**, MC SE **0.0022**.

The shifted-coverage formula is a known-parameter benchmark; it does not claim exact plug-in coverage.

## Verification gates and limitations

Passed **57/57** recorded checks.

Deterministic tolerances and six-standard-error distributional gates are fixed in verify.py. They detect implementation regressions. Six-SE gates are not 95% scientific significance tests. The reported 95% intervals and all seed outcomes are retained whether or not they contain the target.

The simulations are synthetic; agreement does not validate OU assumptions for observed data. Quadrature checks integration and formula implementation; likelihood optimization checks the interior estimator. Neither replaces the derivations.
