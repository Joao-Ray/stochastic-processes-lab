# Supplement: stationary finite-sample calibration and continuous-cell inversion

**Executed methods supplement, 2026-10-10.** This supplements the
[benchmark manuscript](MANUSCRIPT.md). It applies established Monte Carlo
rank-test theory and a standard information bound. It does not establish
priority, an optimal procedure, or journal suitability. Authorship remains
unassigned.

## 1. Question and assumptions

Can inference avoid fitted nuisance-parameter calibration in stationary,
noise-free Gaussian OU records, and can a finite grid support continuous
parameter coverage without interpolation assumptions?

We observe equally spaced values with a stationary random initial state,
unknown equilibrium mean and diffusion scale, and positive drift. Our
calibration is unconditional over that initial state. This differs from the
preceding bootstrap, whose simulated records fix the observed initial value.
Both are evaluated under the same stationary data-generating design; agreement
of their conditional guarantees is not claimed.

## 2. Construction and proof

The [complete proofs](../FINITE_SAMPLE_THEORY.md) establish the following
under ideal Gaussian draws and exact arithmetic:

1. The conditional-profile statistic is invariant to shifting or scaling all
   states. After standardization the stationary law is a canonical AR(1) law,
   depending only on the null slope and record length.
2. Independent null simulation gives the plus-one rank p-value
   $p=(1+\#\{T_j\ge T_0\})/(B+1)$. Exchangeability bounds its rejection
   probability by $\lfloor(B+1)\alpha\rfloor/(B+1)\le\alpha$.
   Conservative tie counting remains valid at the boundary reference zero.
3. For canonical stationary path laws,

$$D(P_\phi\Vert P_\psi)=\frac n2\left[
\log\frac{1-\psi^2}{1-\phi^2}
+\frac{1-\phi^2+(\phi-\psi)^2}{1-\psi^2}-1\right].$$

   On a cell containing the reference slope, the KL supremum occurs at an
   endpoint. Pinsker converts it to a total-variation bound $\delta$.
4. Test each data-independent cell reference at level $\alpha-\delta$,
   retaining its entire cell if not rejected. Cells with $\delta\ge\alpha$
   are retained unconditionally. Under the true slope, the rejection bound is
   $(\alpha-\delta)+\delta\le\alpha$. The union therefore has continuous
   drift-parameter coverage at least $1-\alpha$.

A finite tested grid alone would not have that guarantee. No multiple-testing
correction is required for coverage of a single true parameter: only its
fixed cell matters. Monte Carlo banks must be independent of observed data.
This is an application of [Dufour's Monte Carlo test framework](https://jeanmariedufour.github.io/Dufour_1995_MCT_W.pdf)
and the [standard Pinsker bound](https://web.stanford.edu/class/stats311/lecture-notes.pdf),
not a claim to have invented those tools.

The implemented partition has spacing at most $0.04/\sqrt n$ in
$u=\operatorname{atanh}\phi$, with a final retained cell $[\tanh(5),1)$.
It preserves the small-drift tail. The output is a union, potentially
disconnected, followed by an explicitly conservative enclosing interval.
Its hull always reaches the zero limit; it is not a sharp two-sided interval.

## 3. Numerical design and audit

The [fixed written design](../finite-sample/README.md) and JSON configurations
precede this execution; they are not externally preregistered. The main design
has six cases with 2,000 records each and B=399. Its clean/noisy trajectories
are reused from the original benchmark, with unrestricted slopes checked
against every matched archived record. Reuse enables paired comparisons but
does not add independent observational evidence. The confirmation has 2,000
records in each of two clean cases, new data and bank seeds, and B=999.

Thus **16,000 new inference records** are archived: 12,000 matched-data/main
records (including fixed-start perturbations) and 4,000 independent confirmation
records. Every outer record uses fresh simulation banks. MC standard errors
and Wilson intervals refer to independent outer outcomes within each scenario.
Shared trajectories across main scenarios prevent treating all cases as
independent evidence. The confirmation changes B and seeds together and cannot
isolate either effect.

The measured endpoints are true-point nonrejection and true-cell membership.
Truth is used to measure coverage, never to choose the inference grid or null
calibration law. Full continuous inversions are executed only for three
predetermined replication-zero examples, whose raw series and every cell
decision are saved. This is not a full-set width study.

## 4. Executed findings

### 4.1 Main design

Each row has 2,000 attempted records. Profile, point-rank and cell membership
are available for all records in this run. Coverage percentages are:

| Duration r | Step a | Observation/start | Chi-square profile | Point-rank | TV-buffered cell |
| --- | --- | --- | --- | --- | --- |
| 2 | 0.5 | Clean stationary | 90.35 | 94.55 | 95.50 |
| 2 | 0.1 | Clean stationary | 82.45 | 94.60 | 96.05 |
| 10 | 0.1 | Clean stationary | 90.10 | 94.25 | 95.50 |
| 30 | 0.1 | Clean stationary | 94.50 | 95.40 | 97.05 |
| 10 | 0.1 | Noise SD 0.5 | 10.15 | 15.70 | 19.60 |
| 10 | 0.1 | Fixed X0 at mean + 5 stationary SD | 92.55 | 96.30 | 97.40 |

In the short r=2,a=0.1 design, point-rank minus profile is **12.15 percentage
points**, paired MC SE **0.73 points**; buffered-cell minus profile is
**13.60 points**, paired MC SE **0.77 points**. The reference cell has
TV margin 0.0135 and effective significance level approximately 0.0365.
Results support improved calibration in these clean cases, not uniform
superiority in power or precision.

Old fitted-bootstrap results use different availability populations. At
r=2,a=0.5, bootstrap-profile coverage is 97.87% among only 893 interior fits;
it is not an all-attempt competitor to the 2,000-record rank result. Reports
label that conditional denominator explicitly.

### 4.2 Independent confirmation

| r,a | Profile | Point-rank [Wilson 95%, %] | MC SE, points | Buffered cell [Wilson 95%, %] | MC SE, points |
| --- | --- | --- | --- | --- | --- |
| 2,0.1 | 81.35% | 94.65 [93.58,95.55] | 0.50 | 95.95 [94.99,96.73] | 0.44 |
| 10,0.1 | 91.30% | 95.90 [94.94,96.68] | 0.44 | 96.80 [95.93,97.49] | 0.39 |

The first cell interval overlaps 95%. Neither simulation sample requires its
empirical proportion to exceed the mathematical population lower bound.
This is independent numerical corroboration in two designs, not a proof of
coverage across parameter space. The analytic argument supplies the stated
ideal-model guarantee; the numerical study checks its implementation locally.

### 4.3 Power and full-set informativeness

At r=10,a=0.1, point-rank rejection of false drift 0.5 times truth is
17.25% in the main run and 15.00% in confirmation. Against twice the true
drift it is 21.55% and 22.35%, respectively. These two illustrative points
show limited discrimination at this duration; they are not a full power curve
and do not measure the more conservative cell procedure's power.

Predetermined main examples have:

| r,a | Partition cells | Connected components | Enclosing theta interval |
| --- | --- | --- | --- |
| 2,0.5 | 251 | 1 | [0,infinity] |
| 2,0.1 | 561 | 2 | [0,4.9798] |
| 10,0.1 | 1,251 | 4 | [0,1.9403] |

Intervals denote endpoint closures intersected with positive theta. The first
example is completely uninformative. The last two are unions whose hulls fill
rejected gaps. Finite simulation noise can cause fragmentation; the retained
near-one slope tail can also be isolated. No empirical distribution of width
or component count is inferred from three examples. Figures restrict only
the display range, preserving zero/infinite tails in the saved sets.

### 4.4 Assumption stress

The noise case has buffered membership only 19.60%; valid stationary-null
calibration cannot correct a wrong observation model. The fixed-start case
has high membership in this one run, but violates the stationary-law argument.
It gives neither conditional-on-X0 validity nor a theorem for arbitrary starts.
Both stress cases are retained as diagnostics.

## 5. Reproducibility and scope

[Main report and files](../finite-sample/results/main/report.md) and
[confirmation report and files](../finite-sample/results/confirmation/report.md)
contain all per-record ranks/hits, paired results, MC uncertainty, seeds,
versions, source/configuration hashes and archive SHA256 checksums. Full
example decisions and observations are in the main `examples.json`.
Independent tests check recursion, least-squares likelihood, affine invariance,
rank/tie enumeration, and KL against full Gaussian covariance matrices.
CLI replay checks identical record archives. GitHub Actions repeats both full
runs on Python 3.11, 3.13 and 3.14. Coverage is not a software pass threshold.

Floating-point calculations use stable KL evaluation and a small outward
cushion. This is not certified interval arithmetic for every possible input.
The ideal-model proof is distinguished from pseudorandom numerical execution.

This stage supplies a specified stationary construction, explicit proof and
independent validation. The unanswered research target remains a useful
noise-aware/initial-state-robust procedure with controlled error and strong
published competitors. A focused novelty review, sharper precision analysis,
and stronger power comparisons are still needed before claiming an original
journal contribution.
