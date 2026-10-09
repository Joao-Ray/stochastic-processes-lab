# Finite-window uncertainty in Ornstein–Uhlenbeck drift inference: a reproducible benchmark

**Working manuscript — computational benchmark.** This draft reports executed
simulations and established statistical procedures. It does not claim a novel
estimator, a new calibration theorem, or readiness for submission to a Q1 journal.
Authorship and affiliations have not been assigned.

## Abstract

We examine finite-window inference for the positive drift parameter of a scalar
Ornstein–Uhlenbeck process with unknown equilibrium mean and diffusion scale.
Using exact transitions and paired observation grids, we compare delta-method
Wald, chi-square-profile, basic parametric-bootstrap, and bootstrap-profile
confidence sets across nine clean designs and two measurement-error designs.
The main study contains 22,000 outer records; a separately seeded confirmation
adds 2,000. All failed point fits, unrestricted bootstrap slopes, and unbounded
confidence sets contribute to explicitly defined performance summaries.
With two reversion times and dense sampling, chi-square-profile and bootstrap-
profile conditional coverage are 81.60% and 86.33%, respectively. A short sparse
design produces only 893 admissible point fits out of 2,000; its bootstrap-profile
coverage of 97.87% is accompanied by 867 unbounded upper limits and only 43.70%
all-attempt success-and-coverage. In the central confirmation, bootstrap-profile
coverage is 92.75%, with Wilson interval 91.53%–93.81%. Under substantial
measurement error, the same procedure covers latent drift only 12.10% of the
time. The benchmark distinguishes sampling uncertainty, procedure availability,
and structural misspecification. It provides baselines for further calibration
research rather than evidence of uniformly valid finite-sample inference.

Keywords: Ornstein–Uhlenbeck; profile likelihood; parametric bootstrap;
finite-window inference; measurement error; Monte Carlo uncertainty.

## 1. Introduction and positioning

The OU process is a tractable setting in which an exact transition law does not
automatically give reliable small-sample inference. Its Gaussian transitions
produce a regression representation, but estimated mean reversion can be
strongly affected by elapsed record duration and boundary behavior. Existing
diffusion-bias and bootstrap-correction research includes
[Tang and Chen (2009)](https://doi.org/10.1016/j.jeconom.2008.11.001).
[Yu (2012)](https://doi.org/10.1016/j.jeconom.2012.01.004) studies weak-reversion
bias with a known equilibrium mean; our mean remains unknown.

Parameter uncertainty and sampling design have also been investigated by
[Strey (2019)](https://doi.org/10.1103/PhysRevE.100.062142).
Measurement-noise corrections are addressed in
[Calderon (2013)](https://doi.org/10.1103/PhysRevE.88.012707).
We therefore do not present bootstrap correction, OU interval estimation,
or measurement-error bias as new ideas.

The present objective is narrower: produce an auditable comparison that retains
failures and explains how a high conditional coverage may coexist with frequent
unavailability or uninformative confidence sets. A second objective separates
denser sampling from longer observation using the same latent trajectories.
The protocol and reporting use the simulation-study organization described by
[Morris, White and Crowther (2019)](https://doi.org/10.1002/sim.8086).

## 2. Model and conditional inference

Let

$$dX_t=\theta(\mu-X_t)dt+\sigma dW_t,\quad
\theta>0,\ \sigma>0,\ \mu\in\mathbb R.$$

At spacing $\Delta t$, write

$$X_{k+1}=c+\phi X_k+\varepsilon_k,\qquad
\phi=e^{-\theta\Delta t},\quad c=\mu(1-\phi),\quad
\varepsilon_k\stackrel{\mathrm{iid}}\sim N(0,q),$$

$$q=\frac{\sigma^2}{2\theta}(1-\phi^2).$$

We condition on the observed initial value. For $n$ transitions, the negative
log likelihood is $n\log(2\pi q)/2+S(c,\phi)/(2q)$, where
$S=\sum(x_{k+1}-c-\phi x_k)^2$. The interior solution is ordinary least squares
in intercept and slope, with $\widehat q=S_{\min}/n$. An OU point estimate
exists through this mapping only if $0<\widehat\phi<1$ and $S_{\min}>0$.
Invalid slopes are not converted to successful fits by clipping.

### 2.1 Algebraic profile inversion

For nonconstant lagged values, let
$A=\sum(x_k-\bar x_-)^2$ and let $\widehat\phi$ be the unrestricted slope.
After minimizing the intercept, residual sum of squares satisfies

$$S(\phi)=S_{\min}+A(\phi-\widehat\phi)^2.$$

Profiling $q$ gives a likelihood-ratio statistic
$\Lambda(\phi)=n\log\{S(\phi)/S_b\}$, where $S_b$ is the infimum over
$0<\phi<1$. Its value can be calculated from the nearest boundary when
the unrestricted estimate is outside the domain. This calculates a boundary
supremum; it does not imply finite physical parameters at that boundary.

For any specified cutoff $c_\alpha\ge0$, the profile set is the intersection of
$(0,1)$ with

$$\widehat\phi\pm
\sqrt{\frac{S_b e^{c_\alpha/n}-S_{\min}}A}.$$

**Derivation.** Centering the lagged and following observations eliminates the
intercept. The slope normal equation makes the cross term in the residual
expansion zero. Profiling $q$ reduces the objective to
$n[\log(2\pi S/n)+1]/2$, so the likelihood ratio is a logarithm of residual
sums of squares. Solving $S(\phi)\le S_b e^{c_\alpha/n}$ yields the bounds.
This is standard profile-likelihood algebra, not a newly claimed theorem.

The decreasing transformation $\theta=-\log\phi/\Delta t$ produces the drift
set. A slope lower bound of zero gives an infinite drift upper bound. Endpoint
closures are used for reporting; strictly positive truths have the same membership.

### 2.2 Calibration and correction

The Wald method uses
$\widehat{\mathrm{SE}}(\widehat\theta)
=\sqrt{\widehat q/A}/(\Delta t\widehat\phi)$, with its lower bound truncated
at zero. The chi-square-profile method uses the 95th percentile of
$\chi^2_1$. Both are asymptotic approximations.

Given an interior generating fit, the bootstrap simulates exact fitted-OU
records of the observed length, each starting at the observed first value.
Every unrestricted bootstrap slope is retained. The basic interval reflects
the two slope quantiles about $\widehat\phi$, intersects the physical domain,
and then transforms to drift. The bootstrap-profile method replaces the
chi-square cutoff with the empirical 95th percentile of bootstrap profile
likelihood ratios evaluated at the generating slope. Neither calibration is
asserted to be exact under weak identification.

Separately, we evaluate the one-step correction
$\phi_{bc}=2\widehat\phi-\overline{\widehat\phi^*}$ and map it only when
admissible. This slope correction is distinct from directly correcting theta.
Its performance is compared with the MLE on the common admissible subset.

## 3. Simulation design and analysis

Brownian time and state rescaling reduce the stationary model to drift one,
mean zero, and variance one. We vary dimensionless duration
$r=\theta T\in\{2,10,30\}$ and spacing
$a=\theta\Delta t\in\{0.02,0.1,0.5\}$. At each duration, the three grids
are nested observations of the same fine-grid paths. Initial states follow the
stationary Gaussian law, but their density is excluded from the conditional fit.

The main study uses 2,000 independent outer paths per scenario and $B=399$
bootstrap replicates per admissible fit. At $r=10,a=0.1$, two stress cases
add independent Gaussian observation errors of standard deviation 0.2 or 0.5
to every observation. They share the underlying path and standardized noise.
The fitted procedures intentionally ignore this error and target latent theta.

The central clean scenario is repeated on 2,000 new paths with $B=999$.
This confirmation changes seed and bootstrap size jointly; it is not a controlled
one-factor experiment on $B$. Main and confirmation results are not pooled.

Performance includes conditional coverage, all-attempt success-and-coverage,
availability, empty sets, infinite upper limits, finite mean widths and their
denominators, bias, RMSE, and paired differences. Unavailable methods contribute
zero only to the procedure-level success measure. Coverage MC SE is estimated
as $\sqrt{\widehat p(1-\widehat p)/M}$ using the appropriate denominator;
Wilson intervals accompany coverage fractions. Bias, RMSE, and paired
differences also include MC SE. These are descriptive comparisons, without a
multiplicity-adjusted claim of method superiority.

At nominal 95% coverage, 2,000 outer replications yield MC SE about 0.0049.
The complete [protocol](../PROTOCOL.md), [derivations](../THEORY.md), and
[configuration](../config.json) specify the design and handling of failures.
The protocol was written before these full runs, but was not externally preregistered.

## 4. Results

### 4.1 Duration and calibration

The following clean-design coverage values are conditional on available sets.

| $r$ | $a$ | Wald | Profile | Basic bootstrap | Bootstrap profile |
| --- | --- | --- | --- | --- | --- |
| 2 | 0.02 | 83.58% | 81.60% | 95.56% | 86.33% |
| 10 | 0.02 | 89.80% | 89.40% | 92.45% | 91.80% |
| 30 | 0.02 | 94.20% | 94.15% | 94.05% | 94.85% |
| 30 | 0.1 | 94.90% | 94.50% | 93.55% | 95.20% |

All twelve clean/noisy/confirmation scenarios and uncertainty summaries are
reported in the [main](../results/main/report.md) and
[confirmation](../results/confirmation/report.md) reports. These selected rows
illustrate duration effects, not a selection of a winning method. The basic
bootstrap's favorable short-dense coverage does not extend uniformly across grids.

At $r=10,a=0.1$, the main paired bootstrap-profile minus chi-square-profile
coverage difference is 0.0190 with MC SE 0.0031. At $r=10,a=0.5$ it is
**negative**, -0.0369 with MC SE 0.0047, on their common available subset.
Bootstrap calibration therefore does not uniformly improve coverage.

### 4.2 Availability and unbounded sets

For $r=2,a=0.5$, only 893 of 2,000 unrestricted regressions map to interior
OU estimates. Bootstrap-profile coverage among those records is 97.87%, but
all-attempt success-and-coverage is 43.70%. Its upper endpoint is infinite in
867 of the 893 available sets. The chi-square-profile method is available on
all records, but 1,843 of its 2,000 sets have infinite upper endpoints.

These observations preclude a claim that near-nominal conditional coverage
alone constitutes precise or generally available inference. Finite width
summaries are explicitly conditional on finite nonempty sets.

### 4.3 Point correction and independent confirmation

The correction decreases squared error on the common admissible subsets in
the recorded designs, with paired MSE differences and MC SE provided in the
reports. It can fail: at $r=2,a=0.02$, only 1,412 of 2,000 attempted records
produce an admissible corrected estimate, compared with 1,961 interior MLEs.
This conditional comparison is not a procedure-level guarantee.

In the central confirmation, bootstrap-profile coverage is 92.75%, Wilson
95% interval [91.53%,93.81%]. The paired improvement relative to chi-square-
profile is 0.0280 with MC SE 0.0038, but coverage remains below nominal.
On 1,976 common admissible records, MLE bias/RMSE are 0.5049/0.8471 and
corrected bias/RMSE are 0.0821/0.6887. Thus reducing point bias does not by
itself establish nominal confidence-set coverage.

### 4.4 Observation-model misspecification

With observation error $Y_k=X_k+e_k$, latent stationary variance $V$, and
independent error variance $\eta^2$, the long-record regression pseudo-slope is
$\phi V/(V+\eta^2)$. Consequently

$$\theta_{pseudo}=\theta+\frac{\log(1+\eta^2/V)}{\Delta t}.$$

This is an expanding-domain limit, not the exact mean of a finite-record estimator.
In the noisy central design with $\eta=0.5$, this limit is 3.2314 for true
theta one. The recorded admissible MLE mean is 4.643, still affected by the
short record. Bootstrap-profile and basic-bootstrap cover latent theta only
12.10% and 30.10% of the time. Fitted-model resampling does not correct the
wrong observation model.

![Coverage and noise](../results/main/coverage_and_noise.png)

*Figure 1.* Coverage among available sets, with 1.96 MC SE error bars. The left
panel fixes the finest sampling grid; the right fixes the central duration and
spacing. Availability and width are tabulated separately in the report.

![Point estimation](../results/main/drift_estimation.png)

*Figure 2.* MLE RMSE among interior fits and correction comparisons on their
common admissible subset. Curves across sampling grids use paired latent paths
but may have different admissibility populations; paired MSE comparisons are
reported separately. These are not unconditional procedure risks.

## 5. Discussion and limitations

An exact grid transition removes discretization approximation from this
experiment. It leaves finite-window parameter uncertainty, boundary cases,
and observation-model errors. Dense observations of a short window do not
create a longer elapsed record. Bootstrap-profile calibration can improve a
particular comparison while remaining substantially undercalibrated.

The principal limitations are stationary-start dependence of the simulation
design, limited coverage of alternatives, finite bootstrap quantile accuracy,
and an independent confirmation that changes two factors jointly. There is no
noise-aware state-space fit, finite-sample validity theorem, or comprehensive
comparison against all published calibration methods. The selected literature
review cannot establish novelty. Simulated observations do not constitute
evidence about an actual oceanic, financial, or physical data set.

Further methodological work should state a precise coverage target under weak
reversion and unknown observation noise, establish assumptions and an error
bound or theorem, and benchmark against the strongest relevant published
procedures. The present study supplies reproducible failure cases and baselines
for that work, not a completed original methodological contribution.

## 6. Reproducibility and research materials

The [repository](https://github.com/Joao-Ray/stochastic-processes-lab) contains
the estimator modules, executable design, configuration, tests, main results,
and separate confirmation results. Across both runs, 24,000 outer attempts and
10,282,437 inner bootstrap regressions are evaluated. Every outer record is
saved, including failures; inner samples are reproducible from the recorded
seed scheme. Archives retain inadmissible-slope counts and bootstrap summaries.

JSON results include software versions and SHA256 hashes of configuration,
source, protocol, and compressed records. Figure PDFs preserve vector graphics.
No external observational data or human participants are involved. This draft
and implementation were prepared with Codex assistance; mathematical claims
and manuscript authorship require human verification before submission.

## References

Full machine-readable entries are in [references.bib](references.bib).

1. Tang, C. Y., & Chen, S. X. (2009). Parameter estimation and bias correction for diffusion processes. *Journal of Econometrics*, 149(1), 65–81. [DOI](https://doi.org/10.1016/j.jeconom.2008.11.001).
2. Yu, J. (2012). Bias in the estimation of the mean reversion parameter in continuous time models. *Journal of Econometrics*, 169(1), 114–122. [DOI](https://doi.org/10.1016/j.jeconom.2012.01.004).
3. Calderon, C. P. (2013). Correcting for bias of molecular confinement parameters induced by small time series sample sizes in single-molecule trajectories containing measurement noise. *Physical Review E*, 88, 012707. [DOI](https://doi.org/10.1103/PhysRevE.88.012707).
4. Strey, H. H. (2019). On the estimation of parameters from time traces originating from an Ornstein–Uhlenbeck process. *Physical Review E*, 100, 062142. [DOI](https://doi.org/10.1103/PhysRevE.100.062142).
5. Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation studies to evaluate statistical methods. *Statistics in Medicine*, 38(11), 2074–2102. [DOI](https://doi.org/10.1002/sim.8086).
