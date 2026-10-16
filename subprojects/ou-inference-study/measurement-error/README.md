# Unknown measurement error: inference, identifiability and replication

This stage addresses the failure of noise-free inference when observations
contain independent Gaussian sensor error. It uses established full Gaussian
ML and REML methods; neither likelihood optimization nor their asymptotic
profile cutoff is claimed to give a new finite-sample guarantee.

[Fixed protocol](PROTOCOL.md) · [Derivations and limits](THEORY.md) ·
[Main executed report](results/main/report.md) ·
[Independent confirmation](results/confirmation/report.md) ·
[Untrimmed tails](results/confirmation/tail_report.md) ·
[Numerical amendment](AMENDMENTS.md) ·
[Manuscript supplement](../paper/MEASUREMENT_ERROR_SUPPLEMENT.md) ·
[投稿前研究审查](../SUBMISSION_READINESS_ZH.md)

## Model and questions

Stationary OU latent states have unknown mean, drift and variance. Every site
has k simultaneous readings with independent Gaussian errors of unknown
variance. Is correct noise modeling enough to restore profile coverage? Does
replication help drift inference at the same total scalar-reading budget?
How do weak signal and short elapsed duration limit any possible method?

- Exact joint covariance likelihood, checked against dense matrices and
  sequential Kalman filtering. Noisy observations are not treated as AR(1).
- REML removes the unknown constant mean by residual contrasts.
- Site means and within-site contrasts reduce replicated data without losing
  likelihood information. A banded identity evaluates each objective in O(m).
- Gaussian Fisher information eliminates variance nuisances by a Schur
  complement. Gaussian KL and two-point bounds describe weak-signal limits.

The expected-information benchmark is local and applies to regular unbiased
estimators. It is not a bound on every estimator's finite-sample RMSE. Existing
likelihood/REML/information tools are credited in the derivations and parent
[literature review](../LITERATURE.md); originality is not established.

## Design and outputs

The main design has 14 scenarios, 500 records each. Nine compare 1, 2 and 4
replicates at noise SD 0.5, 1 and 2, with 200 total readings across T=10.
Additional cases vary duration, remove noise or violate assumptions using
colored error and a fixed initial state. Confirmation uses a new seed and
1,000 records each for k=1 and 2 with noise SD 1. Main and confirmation are
separate. Methods are paired within a record; different designs have independent
data, so cross-design differences require independent-sample uncertainty.

Reports retain point-estimate availability, conditional bias/RMSE and MC SE,
paired squared-error differences, all-attempt and available-test coverage,
Wilson intervals, variance estimates and optimizer/white/cap flags. Every
attempt is saved as compressed CSV. The first two records/design/method have
independent differential-evolution audits. Source/configuration/archive hashes
and software versions accompany the records. Internal protocol is not external
preregistration, and optimizer audits are not global certificates.

One initial audit found a weak-identification search discrepancy. The search
was repaired, main data rerun, and a previously viewed confirmation seed was
replaced with unused seed **20261018** before final execution, after a second
diagnostic exposed a low-signal mode missed by both searches. The amendment
retains both diagnostics and separates development results from confirmation.

A white-noise limit has no finite drift estimate; in k=1 data its latent/noise
variance decomposition is unavailable too. Positive sensor noise is required
for replicated fitting; perfectly identical repeats have a singular joint
likelihood and are explicitly rejected. Likelihood search caps are disclosed.
Profile membership is evaluated at truth for coverage; full confidence-set
endpoints and width distributions are not computed.

Numerical conventions include signal fraction <=1e-8 flagged as a white limit
and, for k>1, noise fraction bounded below by 1e-9. A noise-fraction boundary
within 1e-7 is flagged. These are search/diagnostic conventions, not tests of
population identifiability. RMSE MC SE uses a delta approximation and can
itself be unstable when squared errors concentrate in a few records.

## Reproduce

Use the root requirements and Python 3.11 or newer, from the repository root:

```bash
python subprojects/ou-inference-study/measurement-error/run.py
python subprojects/ou-inference-study/measurement-error/run.py --config subprojects/ou-inference-study/measurement-error/confirmation.json --output subprojects/ou-inference-study/measurement-error/results/confirmation
python subprojects/ou-inference-study/measurement-error/tail_analysis.py
python subprojects/ou-inference-study/measurement-error/tail_analysis.py --results subprojects/ou-inference-study/measurement-error/results/confirmation
python subprojects/ou-inference-study/measurement-error/diagnose_extremes.py
```

Use a separate `--output scratch/noisy-test` for trial runs. A reduced
`ci-config.json` supports development; full main and confirmation also run in
GitHub Actions. Coverage is not a software pass threshold. All data here are
synthetic; there is no empirical application or journal acceptance claim.

`tail_analysis.py` is a descriptive follow-up that retains every finite
estimate. `diagnose_extremes.py` reconstructs the worst finite REML error in
each supplied confirmation design, independently checks its objective, and
saves raw data and evaluated profile grids. Those examples are explicitly
post hoc, are not representative-set-width studies, and do not retune inference.
