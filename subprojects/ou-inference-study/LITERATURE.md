# Literature positioning and contribution boundaries

Selected primary research inspected on 2026-10-09. This is a starting bibliography,
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
A proposed method would need a precisely stated claim, a proof or controlled
error bound under explicit assumptions, and independent comparisons against
the strongest relevant published methods. No such claim is established here.
