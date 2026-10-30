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
| Bao, Ullah and Wang, [2014 working paper](https://economics.ucr.edu/repec/ucr/wpaper/201413.pdf), abstract and introduction inspected | Exact finite-sample distribution calculations for an OU mean-reversion estimator, including different initial conditions and drift specifications; interval construction is discussed. | A close comparator for the **noise-free** stage. Match its exact estimator, conditioning and assumptions before claiming a comparison or extension. This review does not establish its publication history. |
| Li, Kenney and Gu, [arXiv:2610.25266](https://arxiv.org/abs/2610.25266), abstract and introduction inspected | OU state-space model with measurement error, likelihood estimation and profile intervals, with microbiome applications and change testing. | Closely overlaps the noisy OU likelihood and profile-inference theme. It is a preprint; its treatment of finite-sample coverage, nuisance parameters, and replicated observations needs a section-by-section comparison before any priority claim. |
| Li, [*On Ornstein-Uhlenbeck State Space Models*](https://dalspaceb.library.dal.ca/server/api/core/bitstreams/f987142c-dea0-4883-a2be-1c1e99e61ce9/content), thesis PDF search-indexed excerpt only | The indexed excerpt mentions sampling designs, replication at time points and a fixed total sample size. | Potentially very close to our replication design. The full PDF timed out; do not claim a precise theorem comparison until the relevant chapter is obtained and read. |
| Chen, Simpson and Ying, [Statistica Sinica paper](https://www3.stat.sinica.edu.tw/statistica/oldpdf/A10n17.pdf), previously inspected and recorded in the [literature map](../LITERATURE.md) | Gaussian OU covariance with independent measurement error and fixed-domain estimation. | Check which covariance combinations are identifiable before proposing uniform precision for separate drift and variance parameters. |

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
2. Verify whether the Bao working paper has a later peer-reviewed version and
   read the relevant derivation and interval-construction sections.
3. Inspect the Li–Kenney–Gu preprint's confidence-set construction and
   simulation designs in full, including boundary and missing-fit handling.
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
