# Independent-design uncertainty for fixed-budget noisy OU inference

Analysis added 2026-10-26, after viewing the final confirmation results.
This retrospective secondary analysis reuses the existing 2,000 final
confirmation records (seed 20261018). It adds no trajectories or likelihood
fits and does not select a replacement seed.

[Executed comparison](results/confirmation/design_comparison.md) ·
[Numerical results](results/confirmation/design_comparison.json) ·
[Original protocol](PROTOCOL.md) · [Archive checks](../../../notes/ARCHIVE_VALIDATION.md)

## Question, units and estimands

At N=200 scalar readings, T=10 and sensor noise SD 1, compare k=1
(200 times) with k=2 (100 times). Repetition changes both the number of
distinct times and noise reduction within a time. It is not a pure
intervention holding the time grid fixed. Other generating parameters
are in `confirmation.json`.

The sampling unit is an entire independently generated synthetic record.
Readings within a record are dependent. Methods fitted to the same record
share data. Different designs use independent simulation streams; matching
their replication indices would create artificial pairs.

For each method and design, define two binary outcomes on **every attempt**:

- A: a finite drift estimate is returned under the archived numerical rules.
- H: a truth-membership test is available and accepts the generating drift.

An unavailable estimate has A=0. An unavailable test has H=0. A and H need
not agree: a boundary fit can permit truth membership without a finite point
estimate. No record is dropped for either endpoint. The targets are
delta(j,E)=P(E=1|k=2)-P(E=1|k=1), E in {A,H}. Four methods and two endpoints
form one fixed family of K=8 contrasts, all reported. H combines numerical
availability with membership; it is not coverage conditional on availability.

## Independent Monte Carlo standard error

For counts xd out of nd records and rates p-hat(d)=xd/nd,

$$
\widehat\delta=\widehat p_2-\widehat p_1,\qquad
\widehat{\mathrm{SE}}_{\rm MC}
=\sqrt{\frac{\widehat p_1(1-\widehat p_1)}{n_1}
       +\frac{\widehat p_2(1-\widehat p_2)}{n_2}}.
$$

Independence across designs removes the covariance term. At an empirical
rate 0 or 1 the plug-in SE can vanish while the population probability is
unknown. No normal interval based solely on that SE is used for the bound.
Normality and equal-variance screens for continuous outcomes are inapplicable
to binary endpoints. Independence comes from the generating design, not a
test on the rows. Archive checks verify complete distinct replication IDs,
binary flags and consistency with likelihood-ratio decisions; they do not
certify the pseudorandom generator.

## Conservative simultaneous bound: derivation

Assume independent records, identically distributed within each design, and
a fixed numerical inference procedure. Apply Hoeffding's bounded-sum
inequality to the weighted independent outcomes:

$$
\widehat\delta=\sum_{i=1}^{n_2}\frac{E_{2i}}{n_2}
                  -\sum_{i=1}^{n_1}\frac{E_{1i}}{n_1}.
$$

The positive terms have range length 1/n2 and the negative terms range
length 1/n1. The sum of squared range lengths is 1/n1+1/n2. Therefore

$$
P\{|\widehat\delta-\delta|\ge t\}
\le 2\exp\!\left[-\frac{2t^2}{1/n_1+1/n_2}\right].
$$

For alpha=.05 choose

$$
t_K=\sqrt{\frac{1/n_1+1/n_2}{2}\log\frac{2K}{\alpha}},\qquad
C_j=[\widehat\delta_j-t_K,\widehat\delta_j+t_K]\cap[-1,1].
$$

Each failure probability is at most alpha/K. A union bound gives
P(all eight true differences belong to their Cj) >= 1-alpha. Dependence
between methods/endpoints within a record is allowed: the union bound does
not require independence between contrasts. Clipping to [-1,1] cannot remove
a true difference. This applies an established inequality; it is not a new
inference theorem for noisy OU parameters.

For n1=n2=1,000, t8=0.07594946, or **7.594946 percentage points**. The bound
is deliberately conservative, especially for rates near 0 or 1. It covers
the fixed family under the stated assumptions, not an outcome-dependent
selection among additional designs, analyses or stopping times. The original
simulation used fixed record counts. This follow-up was not specified before
seeing results and remains retrospective.

The intervals describe Monte Carlo performance probabilities. They do not
repair the OU procedures' asymptotic chi-square cutoff or prove 95%
finite-sample coverage for OU parameters. A higher H can reflect wider sets;
without width results it is not a standalone quality ranking. Inclusion of
zero does not establish equivalence.

## Conditional RMSE and missing estimates

Original RMSE estimates target sqrt(E[(theta-hat-theta)^2 | A=1]) for each
method/design. A depends on both; the finite-estimate populations differ.
The comparison is descriptive and does not establish an unconditional risk
ranking for a procedure that always returns an estimate. Assigning missing
estimates zero error would change the estimand and reward failure. No such
assignment, trimming or winsorization is made. The existing untrimmed tail
audit supplies sensitivity information. This stage does not construct an
interval for the unbounded squared-error contrast.

## Precision planning

For equal counts n in each design, tK=sqrt(log(2K/alpha)/n). Thus

$$
n\ge\left\lceil\frac{\log(2K/\alpha)}{h^2}\right\rceil
$$

is sufficient for an untruncated half width at most h. For K=8, alpha=.05,
h=.02 requires **14,421 attempts per design** using this conservative bound.
It is a sufficient precision budget, not a minimum necessary sample size,
a power calculation or an instruction to execute more simulations. A future
study should fix endpoints, a meaningful effect and a precision/power
criterion before execution. Repeated monitoring needs a sequential guarantee.

## Reproduce and verify

From the repository root, Python >=3.11, no third-party dependencies:

```bash
python scripts/compare_noise_designs.py
python scripts/compare_noise_designs.py --check
python -m unittest discover -s tests -p test_design_comparison.py -v
```

`--output scratch/design-comparison` writes elsewhere. `--check` rejects
absent or stale JSON/Markdown. The command checks the archived evidence
first, then derives counts from records instead of rounded percentages.
Outputs store archive, configuration and analysis-source hashes.
JSON checks allow the archive verifier's numerical tolerance (relative 1e-9,
absolute 1e-11) for cross-platform floating-point differences; counts, types
and hashes must match exactly. Markdown is compared as text.

Tests enumerate every count in four small independent binomial experiments
and check the bound there. Boundary rates, unequal sizes, direction reversal,
invalid inputs, retained failure denominators and stale-report rejection are
also checked. Enumeration checks those implementation cases; the argument
above supplies the general mathematical guarantee.

## References and assistance

Hoeffding, W. (1963). Probability inequalities for sums of bounded random
variables. *Journal of the American Statistical Association*, **58**(301),
13–30. [Publisher DOI](https://doi.org/10.1080/01621459.1963.10500830).
The weighted-sum calculation and family construction above are the project's
application of that inequality.

Analysis and implementation were prepared with AI assistance. The
statistical-analysis skill supplied procedural guidance on estimands,
independence, failures, multiplicity and retrospective labeling. It does
not establish statistical validity or constitute independent peer review.
Software-guidance reference, metadata checked 2026-10-26:
Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026).
*Scientific Agent Skills: A Library of Procedural Knowledge for Research
Agents*. [arXiv:2609.00065](https://doi.org/10.48550/arXiv.2609.00065).
