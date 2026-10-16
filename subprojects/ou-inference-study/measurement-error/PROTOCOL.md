# Fixed design: measurement error, identifiability and equal-budget replication

Written 2026-10-15 before the main/confirmation execution. This is a recorded
internal protocol, not an external preregistration. The code and configurations
are hashed. Any pilot used for software checks is excluded from results.

## Estimand and model

Positive drift theta of a stationary scalar Gaussian OU. Equilibrium mean,
latent stationary variance V and independent Gaussian sensor-noise variance
eta^2 are unknown. At each of m equally spaced sites in [0,T], k measurements
share one latent state but have independent errors. Total budget N=m*k=200
includes both endpoints and all replicates. Larger k means fewer distinct times
and larger time gaps. Cost per scalar sensor reading is held fixed; overheads,
availability of simultaneous replication and timing error are not modeled.

## Designs and methods

Main: 500 records in each of 14 scenarios; seed 20261015. Nine designs cross
k=1,2,4 with noise SD=0.5,1,2 at T=10. T=2 and T=30, k=1, noise SD=0.5 test
duration effects. One T=10 clean design checks the zero-noise boundary.
Two T=10,k=1,SD=0.5 diagnostics violate assumptions: stationary AR(1) sensor
noise with lag coefficient 0.6, and fixed latent X0 at mu+5*sqrt(V).

All use theta=1,mu=0,V=1. Different designs have independent data; methods
within each record are paired. Comparing k does not reuse or pool trajectories.
Latent and sensor randomness use SeedSequence([seed,300,scenario_index,role]),
with role 0 latent and role 1 sensor, separate from earlier stages. Latent paths
use exact transitions, never Euler approximation. For repeats, the conditional
regression baseline receives site averages. The full likelihood receives every
replicate, with sufficient within-site residual sum of squares retained.

Methods: noise-ignorant conditional OLS; stationary full-likelihood ML with noise
fixed zero on site averages (initial-density ablation); exact stationary noisy
Gaussian ML; noisy REML. The last two estimate both covariance shape parameters
by nuisance profiling on a fixed coarse drift scan and three bounded local starts, an additional search
in log innovation/noise coordinates, and an explicit lower-drift-cap profile.
The last two safeguards were added after a numerical audit discrepancy;
[the dated amendment](AMENDMENTS.md) records that repair. No truth enters starting
values or fit selection. Search bounds in theta*dt are [1e-5,20], with the
phi=0 limiting model also checked. Search caps, noise boundaries, white limits,
optimizer flags and likelihood gaps are reported, without converting boundary
limits into finite point estimates. A white-noise limit may yield an uninformative
profile confidence set even when it has no drift point estimate.
For k=1 its latent/noise decomposition is also unidentified; those variance
estimates are unavailable, rather than arbitrary white-limit values.

For ML and REML, membership of the true drift in the profile-LR set uses the
asymptotic chi-square(1) cutoff. It is **not** claimed finite-sample valid. We
compute true-point membership on every record, not full confidence-set endpoints
or widths. A negative likelihood-ratio gap larger than 1e-6 invalidates membership
for that attempt and is recorded. The point estimator is not refit using truth.

## Endpoints and denominators

Every attempt is saved. Report point-estimate availability, conditional bias,
conditional RMSE, and its Monte Carlo SE; matched-subset paired squared-error
differences and MC SE; profile coverage among numerically available tests and
all-attempt success-and-coverage, with binomial MC SE and Wilson intervals.
Neither high coverage nor a lower conditional RMSE is sufficient evidence of
unconditional superiority. Noise estimates and boundary/cap/optimizer counts
are separate outcomes. No arbitrary replacement drift is assigned to failures.

An independent confirmation uses new seed 20261018 and 1,000 records each in
T=10,SD=1,k=1 and k=2. It uses identical algorithms and budget, is not pooled
with the main run, and is not used to tune the method. Earlier seeds 20261016
and 20261017 were inspected during numerical debugging and are not the final
confirmation; the amendment log preserves both diagnostics and repairs.

## Verification and theoretical diagnostics

Independent dense-covariance, Kalman-innovation and residual-contrast references
check the profiled objectives. Expected Fisher information after variance
nuisance elimination examines equal-budget designs; its reciprocal is a local
Cramer-Rao benchmark for regular unbiased estimators, not a universal finite-sample
RMSE bound. A Gaussian KL/two-point argument bounds power and confidence-set
informativeness in a weak-signal submodel. Exact covariance checks validate it.

The first two records of every scenario and both ML/REML fits receive an
independent differential-evolution search (fixed separate seeds), including an
explicit phi=0 comparison. These are numerical optimizer audits, not global
certificates. Discrepancies beyond 1e-4 in profiled objective are retained in
an audit file; software defects must be repaired and all affected results rerun.

Population KL and information diagnostics use declared true parameters; they
are theory illustrations and are never inputs to data fitting. The study does
not assert a universal advantage of replication or claim a new likelihood method.
