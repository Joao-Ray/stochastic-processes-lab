# OU Drift Inference: Bias, Coverage, and Model Misspecification

A reproducible statistical-methods benchmark within Stochastic Processes Lab.
It asks whether nominal 95% confidence sets for the mean-reversion speed remain
calibrated when observation duration is short, sampling is sparse, or the
noise-free observation model is wrong. The equilibrium mean and diffusion
scale are unknown nuisance parameters.

**Status:** completed computational benchmark, stationary finite-sample
calibration extension, and working manuscript. The
methods are established statistical procedures; this project has not established
an original methodological contribution or readiness for a Q1 journal.

[Protocol](PROTOCOL.md) · [Derivations](THEORY.md) · [Main executed report](results/main/report.md) ·
[Independent confirmation](results/confirmation/report.md) ·
[Working manuscript](paper/MANUSCRIPT.md) · [Literature positioning](LITERATURE.md) ·
[中文研究说明](RESEARCH_STATUS_ZH.md)

## Stationary finite-sample extension

The [new stage](finite-sample/README.md) calibrates an affine-invariant profile
statistic using independent stationary-null Monte Carlo ranks, without fitting
nuisance parameters. [Explicit proofs](FINITE_SAMPLE_THEORY.md) combine
exchangeability and a path-KL/Pinsker bound to retain entire parameter cells
with continuous-parameter coverage at least 95% under **stationary, noise-free
Gaussian OU** assumptions. This is an application of existing theory; priority
and novelty have not been established.

The extension archives **12,000 main inference records and 4,000 independent
confirmation records**. Main clean/noisy data are reused and audited against
the preceding study. In the short-record confirmation, chi-square-profile,
point-rank and buffered true-cell coverage are **81.35%, 94.65% and 95.95%**,
respectively. The last Wilson interval is **94.99%–96.73%**. Coverage estimates
include every record, with fresh banks and Monte Carlo uncertainty.

Precision is limited: a sparse example retains all positive drift values;
other examples are disconnected. Full inversions are saved for three fixed
examples only, and every hull reaches zero. Noise-SD-0.5 buffered membership
is only **19.60%**. The guarantee does not extend to noise or fixed initial
states. Two point-null power measurements also show limited discrimination.

[Main report](finite-sample/results/main/report.md) ·
[Independent confirmation](finite-sample/results/confirmation/report.md) ·
[Proofs](FINITE_SAMPLE_THEORY.md) ·
[Manuscript supplement](paper/FINITE_SAMPLE_SUPPLEMENT.md)

## Research design

- Dimensionless duration $r=\theta T\in\{2,10,30\}$ and sampling interval
  $a=\theta\Delta t\in\{0.02,0.1,0.5\}$: nine exact stationary OU scenarios.
- Two observation-error scenarios at $r=10,a=0.1$, with measurement-noise
  SD 0.2 and 0.5 relative to a stationary latent SD of one.
- **2,000 outer paths per scenario**, using 399 parametric-bootstrap paths
  per admissible point fit: **22,000 main attempts**.
- A central-scenario confirmation uses **2,000 new paths and B=999**.
- Methods share each observed series. Coarser grids subsample the same latent
  paths; measurement-error scenarios share latent paths and standardized noise.
- Main and confirmation together execute **10,282,437 bootstrap regressions**.

## Methods

| Method | Calibration | Boundary handling |
| --- | --- | --- |
| Wald | Delta-method conditional regression Hessian | Requires interior OU fit; lower theta bound zero |
| Profile likelihood | Asymptotic chi-square(1) | Boundary likelihood suprema; possibly infinite upper theta |
| Basic bootstrap | Empirical unrestricted AR slope quantiles | All inner slopes retained; intersection with physical domain may be empty |
| Bootstrap profile | Fitted-model likelihood-ratio quantile | All inner slopes retained; requires interior generating fit |

A secondary one-step **slope** bias correction is evaluated separately. It is
not the same as directly correcting theta and is not claimed as a new method.

The reports distinguish coverage among available confidence sets from
success-and-coverage among **all** attempts. They also give finite-width
denominators, unbounded sets, correction failures, paired differences, and
Monte Carlo uncertainty. A high conditional coverage can coexist with a low
probability of producing a useful set.

## Recorded findings

| Scenario | Finding |
| --- | --- |
| Short dense record, $r=2,a=0.02$ | Profile coverage **81.60%**; bootstrap-profile **86.33%** among 1,961 available fits. Basic-bootstrap coverage **95.56%** in that same availability population. This does not establish a universal winner. |
| Short sparse record, $r=2,a=0.5$ | Only **893/2,000** interior fits. Bootstrap-profile conditional coverage **97.87%**, but all-attempt success-and-coverage only **43.70%**. Its upper limit is infinite on **867/893** available records. |
| Long record, $r=30,a=0.1$ | Bootstrap-profile coverage **95.20%**, chi-square-profile **94.50%**. |
| Central confirmation, $r=10,a=0.1$ | Bootstrap-profile coverage **92.75%**, Wilson 95% interval **91.53%–93.81%**; increasing B and changing seed does not remove undercoverage. |
| Noise SD 0.5, $r=10,a=0.1$ | Bootstrap-profile coverage of latent theta **12.10%**; basic-bootstrap **30.10%**. Fitted-model resampling does not repair an incorrect observation model. |

![Coverage and observation error](results/main/coverage_and_noise.png)

Point correction reduces MSE on the common admissible subset in these recorded
scenarios, but can itself be inadmissible. At $r=2,a=0.02$, only **1,412/2,000**
corrected estimates are available. Its conditional RMSE cannot be presented as
the performance of a procedure that always returns an estimate.

## Reproduce

From the repository root, with [requirements.txt](../../requirements.txt) installed:

```bash
python subprojects/ou-inference-study/run.py
python subprojects/ou-inference-study/run.py --config subprojects/ou-inference-study/confirmation.json --output subprojects/ou-inference-study/results/confirmation
```

To try a separate configuration, use `--config` and `--output scratch/my-study`.
A small developer design is available in `ci-config.json`. The full main and
confirmation runs are also executed by GitHub Actions on Python 3.11, 3.13,
and 3.14. Coverage near 95% is not a software pass criterion: poor calibration
is retained as a research result.

Each result folder contains a compressed CSV with **every outer attempt**,
summary JSON, report, PNG figures, and vector PDF figures. Gzip archives can
be read with `gzip.open(..., 'rt')` and `csv.DictReader`. Inner paths and
unrestricted slopes are regenerated from the recorded seed scheme; archives
retain per-record counts, mean slope, calibration cutoff, estimates, and
confidence endpoints. Source, configuration, protocol, and archive SHA256
hashes are included in the JSON.

## Scope

Stationary initialization is part of this design. The likelihood is conditional
on the observed initial state; bootstrap training paths also fix that value.
The confirmation changes both B and seed and therefore cannot isolate their
individual effects. The literature search is a selected positioning review,
not a systematic review or a certification of novelty. The separate extension
adds stationary calibration arguments and a fixed-start diagnostic; it does
not supply noise-aware inference or a general nonstationary-start guarantee.
Originality, sharper precision, and robust calibration remain open work.
