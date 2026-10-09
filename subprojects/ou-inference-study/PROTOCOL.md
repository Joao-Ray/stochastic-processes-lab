# Study protocol: uncertainty in conditional OU drift inference

Status: a computational benchmark protocol, written before the full runs in
this addition. It is **not externally preregistered**, and does not establish
novelty or journal suitability. The configuration and source hashes are saved
with the executed results. No method is selected using the confirmation run.

## Aims

1. Measure coverage, availability, and precision of conditional drift inference
   when only a few mean-reversion times are observed.
2. Separate elapsed observation duration from sampling density using paired,
   nested observations of identical latent paths.
3. Evaluate whether parametric bootstrap calibration or slope bias correction
   improves a chosen performance measure without hiding boundary failures.
4. Expose the failure of these noise-free methods when observations have error.

## Data-generating mechanisms

The latent process is an exact Gaussian OU with stationary initial law.
Canonical values are theta=1, mu=0, sigma=sqrt(2), stationary variance=1.
The design uses dimensionless duration r=theta*T in {2,10,30} and sampling
interval a=theta*dt in {0.02,0.1,0.5}: nine clean scenarios.

At each duration, coarse observations are subsamples of the same fine-grid
paths. Paths are independent across outer replications and duration groups.
The conditional likelihood excludes an initial density even though the
simulation starts at stationarity. The estimator conditions on the initial
observation; it is not treated as a known model parameter. Its stationary
density depends on nuisance parameters and is not claimed to be ancillary.

Two stress scenarios add independent Gaussian measurement errors of SD 0.2
and 0.5 to all observations, including the first, at r=10 and a=0.1. Both
scenarios use the same latent series and the same standardized measurement
noise. Methods intentionally ignore measurement error; their coverage target
is the true latent theta, not the pseudo-true misspecified regression parameter.

## Estimands and methods

Primary estimand: the positive mean-reversion speed theta. Nuisance mu and
sigma are unknown throughout. Secondary targets: raw AR slope and theta bias.

All confidence sets use level 0.95:

- Delta-method Wald interval, truncated at zero, from the conditional
  regression Hessian. Requires an interior OU point fit.
- Profile likelihood interval with chi-square(1) cutoff, including boundary
  suprema and unbounded confidence sets. No clipped point estimate is reported.
- Basic parametric-bootstrap interval in AR slope, transformed after
  intersection with (0,1). Every bootstrap slope is retained before inversion.
- Bootstrap-calibrated profile likelihood: simulate at the interior conditional
  estimate given the observed initial value; use the 95th percentile of the
  bootstrap likelihood ratios to calibrate the observed profile.

Secondary point estimator: theta obtained from corrected slope
phi_bc=2*phi_hat-mean(phi_hat_star), if phi_bc is admissible. This is a one-step
slope correction, **not a reproduction of a published theta correction**, and
not a newly claimed method. Correction failures are recorded without clipping.

The main run uses 2,000 outer replications and B=399 bootstrap paths per
available fit, fixed seed 20261010. A separate central-scenario confirmation
uses 2,000 new paths, B=999, and seed 20261011. This checks stability jointly
across seed and bootstrap size; it does **not isolate** the effect of B.
No invalid bootstrap replicate is omitted or retried. Original fits outside
(0,1) make bootstrap and Wald unavailable; profile intervals remain evaluable.
An empty basic confidence set is available but has coverage zero.

## Performance measures and denominators

Report, for every method and scenario:

- Availability / all attempts, empty-set and infinite-upper-limit fractions.
- Coverage among available sets, with Monte Carlo SE and Wilson interval.
- Success-and-coverage / all attempts (unavailable procedures count as failures);
  this is a procedure-level success measure, not conventional conditional coverage.
- Mean interval width among finite, nonempty sets, and its denominator. Infinite
  widths are separately reported and never silently discarded.
- Bias and RMSE of valid point estimates, Monte Carlo SE, and valid counts.
  MLE/corrected comparisons also use their common admissible subset.
- Paired coverage differences between bootstrap-profile and chi-square-profile
  on common available records; paired MC SE. Joint-success differences are
  reported separately over all attempts.

At true coverage 0.95, 2,000 independent replications give MC SE about 0.0049.
This cannot reliably distinguish coverage differences substantially below one
percentage point. Conditional denominators can be smaller. B=399 uses a noisy
tail quantile; the confirmation run is a limited sensitivity assessment.

## Analysis and reporting rules

Keep all scenarios, all failed outer fits, and all bootstrap inadmissibility
counts. Save one record per outer attempt, interval endpoints, membership,
corrected-estimate availability, seeds, and source/configuration hashes.
Do not choose a winner from a visually favorable scenario. Comparisons are
descriptive with MC uncertainty; there is no family-wise superiority test.
Main and confirmation results are reported separately, with no pooling.

Unit tests and small CI executions check implementation and accounting, not
whether a method reaches 95% coverage. Poor coverage is an experimental result,
not a reason to change a seed or fail a software test.

## Existing literature and scope

Tang and Chen (2009), *Parameter estimation and bias correction for diffusion
processes*, Journal of Econometrics 149, 65–81, already study diffusion-estimator
bias and parametric-bootstrap correction. Their [author-hosted paper](https://www.songxichen.com/Uploads/Files/Publication/Tang-Chen-09-JoE.pdf)
is a background reference; this implementation does not claim their full algorithm.

Yu (2012), *Bias in the estimation of the mean reversion parameter in continuous
time models*, Journal of Econometrics 169, 114–122, studies OU bias with a known
long-run mean, including weak reversion; see the [institutional record](https://smusg.elsevierpure.com/en/publications/bias-in-the-estimation-of-the-mean-reversion-parameter-in-continu/).
Our nuisance mean is unknown, so those formulae cannot be transplanted directly.

Calderon (2013), [arXiv:1304.4196](https://arxiv.org/abs/1304.4196), addresses
OU-related confinement inference with measurement noise. Observation-noise bias
is an established problem, not a novel finding here.

The design/reporting follows the aims, data mechanisms, estimands, methods,
and performance-measures organization of Morris, White and Crowther (2019),
[arXiv:1712.03198](https://arxiv.org/abs/1712.03198), DOI 10.1002/sim.8086.

Before a publishable methodological claim, a broader systematic literature
review, a genuinely new result or procedure, independent validation of it, and
appropriate finite-sample or asymptotic theory remain necessary.
