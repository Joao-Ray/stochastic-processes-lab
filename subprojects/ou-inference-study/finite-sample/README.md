# Finite-sample stationary OU calibration

This stage adds a nuisance-free stationary Monte Carlo rank test and a
conservative continuous-parameter confidence envelope. Its formal arguments
are in [FINITE_SAMPLE_THEORY.md](../FINITE_SAMPLE_THEORY.md).
The parent [benchmark](../README.md) compares approximate fitted-model procedures.

[Executed main results](results/main/report.md) ·
[Independent confirmation](results/confirmation/report.md) ·
[Manuscript supplement](../paper/FINITE_SAMPLE_SUPPLEMENT.md)

## Fixed design before execution

The main design has four clean cases with (theta*T, theta*dt) equal to
(2,0.5), (2,0.1), (10,0.1), and (30,0.1). It also examines the central case
with measurement-noise SD 0.5 and with a fixed initial displacement of five
stationary SDs. These stress cases violate theorem assumptions. Each design
has 2,000 outer records and B=399 fresh null simulations per active test.

The main clean/noisy records reproduce the preceding benchmark's latent
trajectories. Paired comparisons with its bootstrap-profile method audit the
unrestricted slope for every reused record. That reuse is explicitly accounted
for; it is not additional independent observational evidence. A confirmation
uses new data and bank seeds, 2,000 records in each of two clean cases, B=999,
and the same cell resolution. No method is selected using that confirmation.

The key measurements are point-null nonrejection and true-cell membership.
Every outer record has its own independent simulation banks. MC SE and Wilson
intervals use independent outer outcomes. Coverage is not a software pass gate.
False-null rejection is also reported in the central clean design at theta
candidates 0.5 and 2.0 times truth. These are two illustrative power points,
not a full power curve or an adjusted simultaneous testing procedure.

Full continuous inversions are computed only for predetermined replication-zero
examples in the three shorter clean designs. They save every cell decision,
disconnected component, the conservative enclosing interval, and raw observations.
The final untested cell near phi=1 is always included. It is never silently cut
off to make intervals appear shorter. The study does not estimate full-set widths
over all replications.

This written design is not an external preregistration. Configuration and source
hashes accompany the executed results. Existing rank tests and Pinsker bounds
are credited; novelty of the OU-specific envelope has not been established.

## Run

From the repository root, with the root requirements installed:

```bash
python subprojects/ou-inference-study/finite-sample/run.py
python subprojects/ou-inference-study/finite-sample/run.py --config subprojects/ou-inference-study/finite-sample/confirmation.json --output subprojects/ou-inference-study/finite-sample/results/confirmation
```

Output defaults to `results/main`. Use `--output scratch/my-calibration` to
preserve recorded results. Python 3.11 or newer is supported.

## Interpretation

The theoretical coverage lower bound is over both stationary initial states
and simulation banks under a stationary, noise-free Gaussian OU. It does not
guarantee coverage conditional on a particular first observation, or with
unknown measurement error. The envelope can be conservative and disconnected;
its hull adds still more conservatism. Independent exact Gaussian draws and
exact arithmetic enter the proof. The implementation uses pseudorandom draws
and floating-point calculations, with independent numerical benchmarks and
an outward cushion; it is not a universal machine-verified rounding certificate.
