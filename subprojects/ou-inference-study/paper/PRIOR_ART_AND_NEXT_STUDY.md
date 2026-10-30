# Prior-art check and a falsifiable next study

**Date:** 2026-10-30. This is an exploratory, single-reviewer search, not a
systematic review or a certificate of novelty. It was made after the existing
simulation results had been viewed. It contains no new simulation results.

## Question and scope

For a stationary scalar OU process observed with independent Gaussian error,
can a computable confidence set for mean-reversion speed retain stated
finite-sample coverage with unknown nuisance parameters while remaining useful
in a specified short-window, fixed-reading-budget regime?

The question has two separate requirements: coverage and informativeness.
Returning the whole parameter space may satisfy coverage but does not provide
useful precision. This project does **not** currently establish both under
unknown measurement noise. The established noise-free stationary guarantee
must not be extended to this setting without proof.

## Search record and assessed sources

On 2026-10-30, exploratory web searches used the following queries:

- `"Ornstein-Uhlenbeck" "replication" "measurement error" drift estimation sampling design`
- `"On Ornstein-Uhlenbeck State Space Models" Li thesis replication sampling`
- `Ornstein Uhlenbeck measurement noise finite sample confidence interval mean reversion exact inference`

Search ranking is not reproducible database screening; results were selected
for close relevance, with no exhaustive hit count. Sources below are research
papers, a thesis or a preprint; publication status is not inferred from a
search result. See also the earlier [literature map](../LITERATURE.md).

| Source and access | Relevant established material | Relevance and outstanding check |
| --- | --- | --- |
| Bao, Ullah and Wang, [2017 journal article](https://doi.org/10.1080/07474938.2017.1307977), publisher abstract and [2014 working version](https://economics.ucr.edu/repec/ucr/wpaper/201413.pdf) inspected | Exact finite-sample distribution calculations for an OU mean-reversion estimator, including different initial conditions and drift specifications; the working version discusses interval construction. | A close comparator for the **noise-free** stage. Match its exact estimator, conditioning and assumptions before claiming a comparison or extension. The journal article's full text has not been inspected. |
| Li, Kenney and Gu, [arXiv:2610.25266v1](https://arxiv.org/html/2610.25266v1), sections 2.3 and 3.3–3.4 inspected | Noisy OU state-space likelihood, profile intervals, simulation-based coverage and width, plus microbiome applications. | Closely overlaps noisy OU likelihood and profile inference. It is a preprint; replicated-observation handling and numerical failures require further inspection before any priority claim. |
| Li, [*On Ornstein-Uhlenbeck State Space Models*](https://dalspaceb.library.dal.ca/server/api/core/bitstreams/f987142c-dea0-4883-a2be-1c1e99e61ce9/content), thesis PDF search-indexed excerpt only | The indexed excerpt mentions sampling designs, replication at time points and a fixed total sample size. | Potentially very close to our replication design. The full PDF timed out; do not claim a precise theorem comparison until the relevant chapter is obtained and read. |
| Chen, Simpson and Ying, [Statistica Sinica paper](https://www3.stat.sinica.edu.tw/statistica/oldpdf/A10n17.pdf), previously inspected and recorded in the [literature map](../LITERATURE.md) | Gaussian OU covariance with independent measurement error and fixed-domain estimation. | Check which covariance combinations are identifiable before proposing uniform precision for separate drift and variance parameters. |

## Closer comparison: noisy profile likelihood

In sections 3.3–3.4 of Li, Kenney and Gu, the likelihood is evaluated with a
Kalman filter. At fixed drift, the mean and both variance parameters are
optimized; the profile set uses a chi-square(1) cutoff. Their table reports
200 simulated data sets per N, at drift 1, stationary variance 1, measurement
error variance 1 and mean 0. N=100, 500, 1000 and 5000 are reported with
profile coverage of 93%, 96.5%, 95% and 94%, respectively, and median and
interquartile range of interval widths on the log-drift scale. The authors'
section 2.3 separately states a unit time step for its estimator simulations;
do not assume that time-step statement also specifies the coverage table
without further verification.

The preprint states in section 3.3 that the procedure has correct finite-sample
coverage. That is the authors' claim; the chi-square cutoff and the numerical
table alone do not prove uniform finite-sample coverage. A proof under stated
assumptions was not identified in the sections inspected. This is not a claim
that no proof exists elsewhere in the paper or literature.

Our final confirmation uses T=10, a fixed 200-reading budget and two independent
designs (one or two simultaneous readings per time). It records successful
membership at the generating value for ML and REML and reports failures;
it does not compute full noisy confidence-set widths across all records.
The preprint and repository therefore share a closely related likelihood/profile
baseline but have different reported designs and outputs. Their reported
percentages cannot be treated as a head-to-head comparison. A fair comparison
requires matched observation times, initial distribution, nuisance handling,
interval endpoints, optimizer conventions and failure denominators.

For example, 93% coverage from 200 simulations has plug-in Monte Carlo standard
error approximately sqrt(.93*.07/200)=1.8 percentage points. This describes
simulation variability at that design; it is not a confidence guarantee for
the underlying inferential method. No independent replication of the paper's
table has been run here.

## Candidate claim to test, not an established result

Define a compact parameter region with drift and latent variance bounded away
from zero, bounded noise variance, fixed T, and a fully specified observation
grid. The candidate is a finite-sample confidence set with coverage at least
1-alpha uniformly over that region, accompanied by a defensible measure of
set size and an executable algorithm. Whether such a result is new, sharp or
computationally feasible is **unknown**.

To make this falsifiable, a future protocol must specify the exact parameter
region and sampling grid, the nuisance-parameter treatment, numerical error
control, stopping rules, and a nontrivial precision criterion. Before running a
new confirmation, freeze the algorithm and comparison plan. Compare against
the closest published method that can actually be implemented on the same
data, reporting coverage, set width or diameter (including disconnected and
unbounded sets), runtime, and all numerical failures. Use independent seeds
and Monte Carlo uncertainty; do not select a favorable setting after viewing
results. If the method returns the entire search region, count it as a valid
but uninformative set.

The bounded-away-from-zero region is a proposed restricted target, not a claim
of a lower bound for all noise regimes. Analyze approaches to its boundary
separately: weak identification can force intervals to grow. Any extension
outside stationary Gaussian assumptions, fixed grids or bounded nuisance
space requires its own justification.

## Open work before submission

1. Retrieve and extract the thesis chapter, including model, grid, estimands,
   objective and main results; document differences and overlaps explicitly.
2. Read the published Bao article's full derivation and compare its estimator
   and intervals with the earlier working version and this repository.
3. Inspect remaining Li–Kenney–Gu simulation details, including observation
   spacing, boundary handling, missing fits, and reproducibility of width results.
4. Conduct a documented search in appropriate bibliographic databases, with
   inclusion criteria and records of excluded and unavailable reports. This
   exploratory web search cannot establish priority.
5. Formulate and prove a specific result or revise the question if prior work
   already resolves it. Independently review numerical implementation and
   investigate realistic applications only if they support the final claim.

The literature-review skill supplied search and evidence-handling guidance;
AI assistance does not replace domain review. Software-guidance reference:
Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026),
[*Scientific Agent Skills*](https://doi.org/10.48550/arXiv.2609.00065).
