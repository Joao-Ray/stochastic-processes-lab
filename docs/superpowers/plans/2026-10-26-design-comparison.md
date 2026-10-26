# Noisy OU design comparison implementation plan

> Execution: continue the authorized OU inference improvement in this session.
> Use the existing Python project, test-first development and verification before completion.

**Goal:** quantify independent-design Monte Carlo uncertainty for all four
methods in the archived final confirmation, without treating missing estimates
as zero errors or matching unrelated records by replication index.

**Architecture:** a standard-library script verifies the archive, derives eight
binary contrasts (four methods times coverage/finite-point availability), and
writes deterministic JSON and Markdown. A mathematical note proves conservative
simultaneous Hoeffding bounds. Conditional RMSE remains descriptive.

**Design decisions made by the implementer:** use the existing final confirmation
only (k=2 minus k=1); retain every attempt for binary endpoints; use one family
of eight two-sided contrasts at alpha=.05. Main results are already viewed, so
the analysis is explicitly retrospective. A normal approximation would be
narrower but unreliable at availability 0/1; an exact binomial construction
could improve precision at additional implementation cost. Hoeffding provides a
simple independently checkable conservative bound. This is an implementation
choice within the continuing research task, not a user-approved new study.

## Constraints and review focus

- Python >=3.11, no new dependencies for computation or CI.
- Leave pinned inference sources, simulation configurations and archives intact.
- Use independent samples across designs; methods on the same data may depend.
- No claim of unconditional RMSE superiority or equivalence from zero inclusion.
- Reject stale generated reports, altered archive evidence and unsupported designs.
- Check binomial boundary cases, unequal group sizes, invalid counts and absent estimates.

## Task 1: comparison analysis and tests

Files: `scripts/compare_noise_designs.py`, `tests/test_design_comparison.py`.

Interfaces: `contrast(x1,n1,x2,n2,alpha=.05,family=8) -> dict`,
`required_per_design(half_width,alpha=.05,family=8) -> int`,
`build(root=ROOT) -> dict`, `render(result) -> str`.

- [x] Write failing tests for unequal-size independent contrasts, extreme
  counts, invalid inputs, simultaneous widening and exact small-binomial
  coverage enumeration; run the file and confirm missing implementation.
- [x] Implement the bounds, verifier-backed archived-record analysis and CLI
  (`--output`, `--check`); test real confirmed counts and unavailable estimates.
- [x] Test report generation in a temporary folder and stale report rejection.

## Task 2: evidence and publication

Files: `measurement-error/DESIGN_COMPARISON.md`, confirmation comparison outputs,
study/root README links, archive workflow, submission-readiness update.

- [x] Derive the weighted independent-sample Hoeffding bound, family union
  bound and conservative precision budget; distinguish statistical assumptions
  from software checks and MC intervals from OU parameter intervals.
- [x] Generate reports; run all unit tests, archive verification and report check.
- [x] Review the diff, including the fixed family, independent sample units,
  failure denominators, pinned source files and generated outputs.

Verification: the full suite collected before the final roundoff test passed
155 tests and 11 subtests; the final comparison file passed all 10 tests after
that addition. Archive verification passed all 9,000 records; report check
and local Markdown-link checks passed. A separate author self-review found
no remaining material issue; this is not independent scientific review.
JSON comparisons allow the archive checker's floating tolerance to avoid
false drift across libm versions, while integer counts and hashes remain exact.

Integration gate: commit and publish the verified tree, then verify the remote
head/tree and the new archive CI. A pending full CI must be reported as pending.
