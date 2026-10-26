# Independent-design comparison: final confirmation

Retrospective secondary analysis of the existing 2,000 synthetic records, seed 20261018. No new inference attempts. Both designs spend 200 scalar readings over T=10 with noise SD 1; k=1 uses 200 times, k=2 uses 100 times.

[Derivation and assumptions](../../DESIGN_COMPARISON.md) · [Original report](report.md) · [Untrimmed tails](tail_report.md) · [Machine-readable analysis](design_comparison.json)

## All-attempt binary contrasts

Direction is **k=2 minus k=1**. Rates include all 1,000 attempts/design. An unavailable truth-membership test contributes zero to successful truth membership. Finite-point availability is a separate endpoint; its denominator is also 1,000. Designs use independent data: equal replication indices do not create pairs.

Intervals jointly cover the eight fixed probability differences with probability at least 95% under independent simulation records (Hoeffding plus union bound). They quantify Monte Carlo uncertainty at these designs. They do not establish 95% coverage for OU parameter intervals. MC SE is the plug-in independent-sample standard error; at rates 0/1 it can be zero although uncertainty remains.

| Method | Endpoint | k=1 successes/n | k=2 successes/n | Difference (pp) | MC SE (pp) | Simultaneous 95% bound (pp) |
| --- | --- | --- | --- | --- | --- | --- |
| Noise-ignorant OLS | Successful truth membership | 0/1000 | 10/1000 | +1.00 | 0.31 | [-6.59, +8.59] |
| Noise-ignorant OLS | Finite drift available | 1000/1000 | 1000/1000 | +0.00 | 0.00 | [-7.59, +7.59] |
| Noise-ignorant stationary | Successful truth membership | 0/1000 | 8/1000 | +0.80 | 0.28 | [-6.79, +8.39] |
| Noise-ignorant stationary | Finite drift available | 994/1000 | 992/1000 | -0.20 | 0.37 | [-7.79, +7.39] |
| Noisy ML | Successful truth membership | 930/1000 | 914/1000 | -1.60 | 1.20 | [-9.19, +5.99] |
| Noisy ML | Finite drift available | 998/1000 | 999/1000 | +0.10 | 0.17 | [-7.49, +7.69] |
| Noisy REML | Successful truth membership | 957/1000 | 946/1000 | -1.10 | 0.96 | [-8.69, +6.49] |
| Noisy REML | Finite drift available | 973/1000 | 983/1000 | +1.00 | 0.66 | [-6.59, +8.59] |

0 of the eight conservative intervals exclude zero. Zero inclusion does not prove equivalence or absence of an effect. Higher membership alone also does not imply better calibration: interval width and departure from the nominal target matter.

## Conditional RMSE: descriptive comparison

| Method | k=1 finite n | k=1 RMSE | k=2 finite n | k=2 RMSE |
| --- | --- | --- | --- | --- |
| Noise-ignorant OLS | 1000 | 19.4352 | 1000 | 6.5930 |
| Noise-ignorant stationary | 994 | 19.4257 | 992 | 6.5843 |
| Noisy ML | 998 | 1.7563 | 999 | 1.5114 |
| Noisy REML | 973 | 1.5394 | 983 | 1.2986 |

These RMSE values condition on method- and design-specific finite-estimate events. Missing boundary or optimizer estimates are retained as failures in availability, never assigned zero squared error. The values cannot establish an unconditional risk ranking. Squared-error tails are untrimmed; the [tail audit](tail_report.md) describes their influence. No cross-design RMSE significance claim or paired-record test is made.

## Conservative Monte Carlo precision budget

For the same eight-contrast family and equal sample sizes, the sufficient per-design count for an untruncated half width h is ceil(log(320)/h²). This distribution-free bound is conservative; these are proposed counts, not executed studies, minimum necessary counts or power calculations.

| Target half width (pp) | Sufficient attempts per design |
| --- | --- |
| 5 | 2,308 |
| 2 | 14,421 |
| 1 | 57,684 |

This analysis was specified after inspecting the earlier results and is not a newly preregistered confirmation. Its bounds cover this fixed eight-item family; they do not cover data-dependent selection among all past/future designs, repeated monitoring or all claims in the repository. The original inference algorithms, archived records and seeds are unchanged.
