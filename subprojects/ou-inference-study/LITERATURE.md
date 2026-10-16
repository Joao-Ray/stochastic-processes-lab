# Literature positioning and contribution boundaries

Selected primary research inspected on 2026-10-09. This is a starting bibliography,
Additional calibration sources inspected on 2026-10-10. This is
not a systematic search, and it cannot establish that an unobserved competitor
or theorem does not exist. No journal quartile is asserted here.

| Source | What is already established | Relationship to this repository |
| --- | --- | --- |
| Tang & Chen (2009), [DOI 10.1016/j.jeconom.2008.11.001](https://doi.org/10.1016/j.jeconom.2008.11.001), [author-hosted paper](https://www.songxichen.com/Uploads/Files/Publication/Tang-Chen-09-JoE.pdf) | Diffusion parameter bias/variance analysis and parametric-bootstrap bias correction | Bootstrap correction is prior work. The implemented slope correction is explicitly described; it is not claimed to reproduce their full procedure. |
| Yu (2012), [DOI 10.1016/j.jeconom.2012.01.004](https://doi.org/10.1016/j.jeconom.2012.01.004), [institutional abstract](https://smusg.elsevierpure.com/en/publications/bias-in-the-estimation-of-the-mean-reversion-parameter-in-continu/) | OU mean-reversion bias, including weak reversion, with known long-run mean | This benchmark estimates the mean, so the paper's known-mean formulas are not inserted without derivation. |
| Calderon (2013), [arXiv:1304.4196](https://arxiv.org/abs/1304.4196), [DOI 10.1103/PhysRevE.88.012707](https://doi.org/10.1103/PhysRevE.88.012707) | OU-related confinement inference and finite-sample correction with measurement noise | Noise-induced bias is established. This repository demonstrates the limitation of noise-free intervals; it does not supply a new noise-aware estimator. |
| Strey (2019), [arXiv:1805.05977](https://arxiv.org/abs/1805.05977), [DOI 10.1103/PhysRevE.100.062142](https://doi.org/10.1103/PhysRevE.100.062142) | Likelihood/Bayesian treatment of OU parameter uncertainty and sampling design | Accurate OU intervals and sampling optimization are existing research. His fixed-number-of-samples design differs from paired grids at fixed elapsed duration here. |
| Morris, White & Crowther (2019), [arXiv:1712.03198](https://arxiv.org/abs/1712.03198), [DOI 10.1002/sim.8086](https://doi.org/10.1002/sim.8086) | Structured simulation design and reporting, including Monte Carlo uncertainty | Motivates the protocol, denominators, and MC SEs; reporting rigor is not itself proof of a new statistical method. |
| Liu, Feng & Xiao (2026), [DOI 10.1007/s10463-026-00993-w](https://doi.org/10.1007/s10463-026-00993-w) | Fixed-domain covariance estimation for anisotropic OU fields in dimension at least two | A recent adjacent result. Its field model and estimands differ from scalar OU conditional drift inference; its theorems are not imported here. |

## What this addition provides

### Stationary calibration extension

[Dufour (2006)](https://doi.org/10.1016/j.jeconom.2005.06.007),
[author-hosted working version](https://jeanmariedufour.github.io/Dufour_1995_MCT_W.pdf),
provides the general Monte Carlo test framework, including finite-sample rank
calibration and nuisance-parameter issues. The plus-one p-value and
exchangeability principle in this extension are established tools.
[Stanford information-theory notes](https://web.stanford.edu/class/stats311/lecture-notes.pdf)
give the standard Pinsker inequality with natural-log KL divergence.

The [OU specialization](FINITE_SAMPLE_THEORY.md) uses an affine-invariant
statistic to remove stationary location/scale nuisance parameters, derives
the canonical path KL, and buffers cell tests by a total-variation bound.
The derivation and executed verification are explicit. A selected search
cannot establish priority for this envelope. Its noise-free stationary
assumptions and limited precision leave the broader research target open.

### Preceding benchmark

- A fully specified, reproducible comparison of four uncertainty procedures.
- Algebraic profile inversion with explicit empty/unbounded/boundary handling.
- Paired observation designs, fit-failure accounting, and MC uncertainty.
- A separate confirmation run that retains the observed undercoverage.
- A demonstrable distinction between finite-sample bias and observation-model
  misspecification, with both theory and numerical evidence.

These are useful benchmark and software contributions. None is currently
claimed as an original statistical theorem, a novel estimator, or a result
that by itself meets the contribution standard of a particular journal.

## Candidate research question to resolve next

Can a calibration procedure maintain specified drift-interval coverage over
short observation windows **and** accommodate unknown observation noise, while
retaining boundary cases and quantifying its finite-sample error?

Answering this requires a focused search of simulation-based test inversion,
state-space OU inference, weak-identification theory, and bootstrap calibration.
The present chi-square and plug-in bootstrap procedures provide baselines.
A broader proposed method would need a precisely stated claim, a proof or controlled
error bound under explicit assumptions, and independent comparisons against
the strongest relevant published methods. The stationary envelope addresses
a restricted noise-free case; the broader noise-aware claim is not established.

## Measurement-error extension: additional primary sources (2026-10-15)

| Source | Inspected material | Consequence for contribution claims |
| --- | --- | --- |
| Chen, Simpson & Ying (2000), [journal PDF](https://www3.stat.sinica.edu.tw/statistica/oldpdf/A10n17.pdf), Statistica Sinica 10, 141–156 | Full paper, especially its Gaussian covariance law and infill estimands | OU plus independent white noise is established. Its fixed-domain identifiable product is V*theta; its results do not imply consistency of unknown drift and V separately. The PDF corrects an author-name error on the journal HTML abstract page. |
| Holý & Tomanová, [arXiv:1811.09312](https://arxiv.org/html/1811.09312) (inspected v4, 2025) | Full available text, covariance moments and ARMA/likelihood discussion | Noise correction and an ARMA representation are existing work. This repository uses the full joint covariance density, verified by filtering, and does not claim to reproduce every estimator/experiment in that paper. |
| Patterson & Thompson (1971), [DOI](https://doi.org/10.1093/biomet/58.3.545) | Primary journal bibliographic record; original article text not accessed here | REML is credited as an established variance-component approach. The implemented Gaussian derivation is supplied and checked independently; it is not a new estimator principle. |
| Li, Kenney & Gu (2026), [arXiv:2605.01262](https://arxiv.org/html/2605.01262) | Full available preprint text, stationary Kalman likelihood and identifiability sections | Noisy OU state-space models, multivariate identifiability and SST applications already exist. Our scalar fixed-budget benchmark neither reproduces nor supersedes their factor-model study. |
| Li (2025), [Dalhousie thesis](https://dalspaceb.library.dal.ca/server/api/core/bitstreams/f987142c-dea0-4883-a2be-1c1e99e61ce9/content), *On Ornstein-Uhlenbeck State Space Models* | Institutional metadata and search-indexed text identify sampling/replication work; full PDF retrieval failed in this session | This is a close priority-check lead. Replication under noisy OU must not be presented as a newly discovered topic. Full chapter comparison remains required. |

The [new derivations](measurement-error/THEORY.md) use Gaussian likelihood,
within-site contrasts, a Fisher-information Schur complement and KL/TV
two-point arguments. Those tools are standard. A selected literature search
still cannot certify a novel theorem, algorithm or sampling design. The
[submission-readiness review](SUBMISSION_READINESS_ZH.md) records remaining
contribution and external-review gaps.
