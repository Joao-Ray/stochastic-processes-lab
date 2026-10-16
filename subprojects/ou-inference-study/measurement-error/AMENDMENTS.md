# Numerical amendment log

## 2026-10-15: search near the REML lower-drift cap

The initial main run's independent optimizer audit found one objective
discrepancy exceeding the fixed 1e-4 tolerance. At T=10,k=2,noise SD=1,
replication 0, REML's local objective was 291.4423641163105; differential
evolution found 291.43653156904224 near the lower-drift search cap, a gap of
0.005832547268255439. The [initial audit](development/initial_optimizer_audit.json)
is retained. Initial results were software-development output and are excluded
from the final scientific tables and record count.

Repair: retain the existing three-start scan, add a local search in
(log(theta*dt), log(eta^2/[V*(1-phi^2)])) coordinates, and profile the
lower-drift cap explicitly. These steps use the observed data only. Cap
solutions remain unavailable finite drift estimates. No profile truth value
is used to initialize or repair the point fit. Likelihood-ratio gaps and
nonconvergence remain explicit per-record failure outcomes, not substituted
successes. The change improves a numerical implementation of the same ML/REML
objectives, and is not claimed as a new statistical method.

The seed-20261016 confirmation had already been viewed during debugging.
It is therefore excluded as confirmatory evidence. The final confirmation
configuration was changed to unused seed **20261017** before its execution,
keeping the two designs, sample count and budget fixed. Main seed 20261015
is retained as development/main data. Both final main and confirmation are
rerun with the repaired source, with separate source/configuration hashes.

The final study reports 7,000 main attempts plus 2,000 final confirmation
attempts. Discarded pilots, the initial run and the superseded confirmation
are not counted as additional independent evidence. This internal amendment
log is not an external preregistration.

## 2026-10-15: low-signal competing likelihood mode

A post hoc tail-profile check of seed-20261017 confirmation replication 416,
k=1, revealed a lower REML profile objective at a moderate drift than the
recorded extreme estimate 82.1048. Both the original local fit and the separate
differential-evolution search had selected the same inferior mode. Their
agreement was therefore insufficient evidence of global correctness. The
[initial profile diagnostic](development/initial_extreme_profiles.json) retains
the negative profile-LR values and raw observations. That confirmation is
also software-development evidence and is not used in final scientific tables.

Cause and repair: a coarse joint scan can miss a narrow low-signal region in
the variance fraction. Each fixed coarse drift now profiles the nuisance
fraction before selecting local-search starts. The same likelihood objectives,
three local starts, transformed search, caps and failure conventions remain.
This change is numerical, not a statistical calibration proposal. Regression
tests retain both discovered examples.

The final untouched confirmation seed is **20261018**. It supersedes both
20261016 and 20261017, which were viewed during debugging. Final main and
confirmation are rerun with this source. Post hoc diagnostics may inspect the
new records but do not change scientific fits, and any remaining failures are
explicitly reported. Only the final 7,000 + 2,000 attempts enter scientific tables.
