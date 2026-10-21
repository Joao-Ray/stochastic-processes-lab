# Independent validation of the noisy-OU result archives

Run from any directory with Python 3.11 or later; no third-party packages are
required. From the repository root:

```bash
python scripts/verify_noise_results.py
python -m unittest discover -s tests -p test_result_integrity.py -v
```

The first command checks the committed 7,000-record main archive and the
2,000-record confirmation archive. It does not regenerate or overwrite either
archive. `--output /tmp/noisy-validation.json` optionally saves a machine-readable
receipt containing the verifier, configuration and archive SHA-256 values.
The committed [receipt](noisy_archive_validation.json) records the current
validation; rerunning the command checks whether it still holds.

## What is checked

The verifier is separate from the study runner. It never imports the inference
modules, the summary generator or the tail-analysis script. Calculations use
Python's `math` and `statistics` standard-library modules. Tolerance is 1e-9
relative or 1e-11 absolute for recomputed floating summaries; integer counts,
flags and hashes must match exactly.

| Evidence | Checks |
| --- | --- |
| Provenance | Required source-file list and hashes, configuration content and hash, compressed archive hash |
| Experiment population | Every configured scenario, the complete replication index range, duplicate/missing IDs, malformed CSV and nonfinite estimates |
| Fit status | Positive available drift, drift/slope relation, variance decomposition, white-limit and optimizer status, membership availability and recorded LR cutoff |
| Reported estimation error | Finite-estimate denominators, bias, RMSE, bias MC SE and delta-method RMSE MC SE |
| Reported membership | Available-test and all-attempt denominators, proportions, binomial MC SE and Wilson intervals |
| Method contrasts | Common finite-estimate population, paired squared-error difference and its MC SE |
| Diagnostics | White limits, search caps, noise boundaries, nonconvergence, LR gaps, exceptions, noise-estimate counts and means |
| Tails | Median, linearly interpolated 99th percentile, maximum absolute error, top-five squared-error share, record identity and rank |
| Planned optimizer checks | Every planned successful fit has an audit; fixed audit seeds, local objective, alternative objective and discrepancy flags agree |
| Selected extremes | Post hoc label, configuration/source/archive hashes, worst finite REML record, raw-data dimensions, profile-grid structure and recorded numerical gaps |

Missing membership counts as failure for the all-attempt proportion. A missing
finite point estimate does not automatically imply missing membership. The
checker preserves these distinct populations. It does not require measured
coverage to exceed 95%, or remove records with large errors or failed fits.
An unresolved planned optimizer discrepancy above 1e-4 or selected-profile gap
below -1e-6 fails numerical validation; these are software checks, not an
acceptance threshold for scientific performance.

The corruption tests change a copied bundle: they alter RMSE, coverage
intervals, paired denominators, tail shares, audit evidence and provenance;
remove or duplicate records; insert NaN; and change hit/availability flags.
For CSV changes they also refresh all archive checksums. This checks that the
verifier detects inconsistent results even when the checksum is up to date.
The tests leave the research archives untouched.

## What a passing check does not establish

- It does not establish that synthetic records came from the claimed population.
  Reconstructing the exact seeded observations and checking numerical fits is
  part of experiment reproduction and the selected-record diagnostics.
- It does not independently recompute Fisher matrices, KL bounds or every
  likelihood. Separate dense-matrix, filtering, contrast and analytic tests
  cover those calculations in the full test suite.
- For selected extremes it checks the saved grid and its connection to record
  IDs; it does not recompute the raw observations or profile likelihood grid.
- It does not certify global likelihood optima, valid model assumptions,
  uniform confidence-set coverage, a novel method or publication readiness.
- Source hashes detect drift from this recorded run. They are not a signature,
  external timestamp or proof of authenticity. If the scientific source changes,
  the old receipt will fail against the changed source until the affected study
  is reproduced. Do not merely replace hashes to make it pass.

The [Archive integrity workflow](../.github/workflows/archive-integrity.yml)
runs these checks on Python 3.11, 3.13 and 3.14 without scientific-package
installation. The existing [Tests workflow](../.github/workflows/tests.yml)
continues to execute the full numerical experiments. A fast consistency check
and a new numerical execution answer different questions.

## Recorded verification of the previous scientific commit

The [full numerical workflow for commit 3ffe120](https://github.com/Joao-Ray/stochastic-processes-lab/actions/runs/38842613851)
completed successfully on all three Python versions. That is verification of
that scientific source revision, not a statement about future commits or a
statistical proof. This archive checker adds a separately runnable check of
the published records and reported summaries.
