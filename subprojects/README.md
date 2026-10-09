# Research Subprojects

These studies have their own research question, configuration, executable entry
point, numerical record, figures, and report. They reuse the parent repository's
mathematical modules and Python environment.

| Subproject | Research question | Entry point | Executed results |
| --- | --- | --- | --- |
| [OU Calibration and Forecasting](ou-calibration/README.md) | How do observation design, estimated parameters, and a future mean shift affect recovery and prediction? | `python subprojects/ou-calibration/run.py` | [Report](ou-calibration/results/report.md) |
| [OU Drift Inference Study](ou-inference-study/README.md) | How do finite windows, boundary failures, and observation noise affect drift confidence-set coverage? | `python subprojects/ou-inference-study/run.py` | [Main](ou-inference-study/results/main/report.md), [confirmation](ou-inference-study/results/confirmation/report.md) |

The first subproject uses synthetic data with known truth so estimation and
forecast errors can be measured. Future empirical applications should report
their data source, preprocessing, and validation separately.

The OU study also includes [theoretical derivations](ou-calibration/THEORY.md)
and a separate [executed numerical verification](ou-calibration/results/verification/report.md).
Run `python subprojects/ou-calibration/verify.py` to reproduce its quadrature,
optimization, moment, MSE, and mean-shift comparisons.

The inference study adds a [fixed design protocol](ou-inference-study/PROTOCOL.md),
[uncertainty derivations](ou-inference-study/THEORY.md), and
[working manuscript](ou-inference-study/paper/MANUSCRIPT.md). Its methods are
established baselines; a novel statistical contribution is not yet claimed.
