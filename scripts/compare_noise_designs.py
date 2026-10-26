"""Retrospective independent-design comparisons of final noisy-OU confirmation.

Only archived data are used. No new simulations or likelihood fits are run.
The simultaneous bound concerns Monte Carlo performance probabilities, not
confidence sets for OU parameters. See measurement-error/DESIGN_COMPARISON.md.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import verify_noise_results as verifier

METHODS = verifier.METHODS
ENDPOINTS = ("success_and_coverage", "finite_point_available")
FAMILY = len(METHODS) * len(ENDPOINTS)
ALPHA = .05
DEFAULT_OUTPUT = ROOT / verifier.STUDY / "results/confirmation"


def settings(alpha, family):
    if not isinstance(alpha, (int, float)) or isinstance(alpha, bool) or not 0 < alpha < 1:
        raise ValueError("alpha must be in (0,1)")
    if type(family) is not int or family < 1:
        raise ValueError("family must be a positive integer")


def contrast(x1, n1, x2, n2, alpha=ALPHA, family=FAMILY):
    """Return p2-p1 with independent MC SE and simultaneous Hoeffding bound."""
    settings(alpha, family)
    for hits, size in ((x1, n1), (x2, n2)):
        if type(size) is not int or size < 1 or type(hits) is not int or not 0 <= hits <= size:
            raise ValueError("counts must be integers with 0 <= hits <= n and n > 0")
    p1, p2 = x1/n1, x2/n2
    difference = p2-p1
    width = math.sqrt(.5 * (1/n1 + 1/n2) * math.log(2*family/alpha))
    return dict(first=dict(hits=x1, n=n1, rate=p1), second=dict(hits=x2, n=n2, rate=p2),
                difference=difference, mc_se=math.sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2),
                half_width=width,
                simultaneous_interval=[max(-1., difference-width), min(1., difference+width)])


def required_per_design(half_width, alpha=ALPHA, family=FAMILY):
    """Conservative equal-size precision budget; not a power calculation."""
    settings(alpha, family)
    if (not isinstance(half_width, (int, float)) or isinstance(half_width, bool)
            or not 0 < half_width <= 1):
        raise ValueError("half width must be in (0,1]")
    return math.ceil(math.log(2*family/alpha) / half_width**2)


def build(root=ROOT):
    root = Path(root)
    base = root / verifier.STUDY
    folder = base / "results/confirmation"
    receipt = verifier.verify_bundle(folder, base / "confirmation.json", root, require_extremes=True)
    cfg = verifier.read_json(base / "confirmation.json")
    designs = [dict(duration=10., repeats=k, noise_sd=1., stress="none") for k in (1, 2)]
    # This report has a deliberately fixed scope, including the viewed seed.
    verifier.compare(cfg, dict(seed=20261018, replications=1000, budget=200, theta=1., mu=0.,
                               stationary_variance=1., level=.95, audit_records_per_scenario=2,
                               scenarios=designs), "unsupported confirmation configuration")
    ids = [f"T=10,k={k},eta=1,stress=none" for k in (1, 2)]
    with gzip.open(folder / "replications.csv.gz", "rt", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    blocks = [[row for row in rows if row["scenario"] == name] for name in ids]
    contrasts, risks = [], []
    for method in METHODS:
        for endpoint in ENDPOINTS:
            counts = []
            for block in blocks:
                if endpoint == "success_and_coverage":
                    hits = sum(verifier.boolean(row, method + "_test_available")
                               and verifier.boolean(row, method + "_hit") for row in block)
                else:
                    hits = sum(row[method + "_theta"] != "" for row in block)
                counts.append(hits)
            contrasts.append(dict(method=method, endpoint=endpoint,
                                  **contrast(counts[0], len(blocks[0]), counts[1], len(blocks[1]))))
        points = []
        for block in blocks:
            estimates = [verifier.number(row, method + "_theta") for row in block
                         if row[method + "_theta"] != ""]
            points.append(verifier.point_metrics(estimates, cfg["theta"]))
        risks.append(dict(method=method, first=points[0], second=points[1]))
    sources = [Path(__file__).resolve().relative_to(ROOT), Path("scripts/verify_noise_results.py"),
               verifier.STUDY / "DESIGN_COMPARISON.md"]
    return dict(analysis="retrospective; archived final confirmation only; no new simulation",
                direction="k=2 minus k=1", seed=cfg["seed"], budget=cfg["budget"],
                alpha=ALPHA, family_size=FAMILY, designs=ids, archive_validation=receipt,
                source_sha256={path.as_posix(): verifier.digest(root / path) for path in sources},
                contrasts=contrasts, conditional_rmse=risks,
                precision_budget=[dict(half_width=h, replications_per_design=required_per_design(h))
                                  for h in (.05, .02, .01)])


def render(result):
    labels = {"conditional": "Noise-ignorant OLS", "stationary_ignorant": "Noise-ignorant stationary",
              "ml": "Noisy ML", "reml": "Noisy REML"}
    endpoints = {"success_and_coverage": "Successful truth membership",
                 "finite_point_available": "Finite drift available"}
    lines = ["# Independent-design comparison: final confirmation", "",
             "Retrospective secondary analysis of the existing 2,000 synthetic records, "
             "seed 20261018. No new inference attempts. Both designs spend 200 scalar readings "
             "over T=10 with noise SD 1; k=1 uses 200 times, k=2 uses 100 times.", "",
             "[Derivation and assumptions](../../DESIGN_COMPARISON.md) · "
             "[Original report](report.md) · [Untrimmed tails](tail_report.md) · "
             "[Machine-readable analysis](design_comparison.json)", "",
             "## All-attempt binary contrasts", "",
             "Direction is **k=2 minus k=1**. Rates include all 1,000 attempts/design. "
             "An unavailable truth-membership test contributes zero to successful truth membership. "
             "Finite-point availability is a separate endpoint; its denominator is also 1,000. "
             "Designs use independent data: equal replication indices do not create pairs.", "",
             "Intervals jointly cover the eight fixed probability differences with probability "
             "at least 95% under independent simulation records (Hoeffding plus union bound). "
             "They quantify Monte Carlo uncertainty at these designs. They do not establish "
             "95% coverage for OU parameter intervals. MC SE is the plug-in independent-sample "
             "standard error; at rates 0/1 it can be zero although uncertainty remains.", "",
             "| Method | Endpoint | k=1 successes/n | k=2 successes/n | Difference (pp) | MC SE (pp) | Simultaneous 95% bound (pp) |",
             "| --- | --- | --- | --- | --- | --- | --- |"]
    for row in result["contrasts"]:
        a, b = row["first"], row["second"]
        lo, hi = row["simultaneous_interval"]
        lines.append(f"| {labels[row['method']]} | {endpoints[row['endpoint']]} | "
                     f"{a['hits']}/{a['n']} | {b['hits']}/{b['n']} | "
                     f"{100*row['difference']:+.2f} | {100*row['mc_se']:.2f} | [{100*lo:+.2f}, {100*hi:+.2f}] |")
    resolved = sum(not lo <= 0 <= hi for lo, hi in (row["simultaneous_interval"] for row in result["contrasts"]))
    lines += ["", f"{resolved} of the eight conservative intervals exclude zero. "
              "Zero inclusion does not prove equivalence or absence of an effect. Higher "
              "membership alone also does not imply better calibration: interval width and "
              "departure from the nominal target matter.", "",
              "## Conditional RMSE: descriptive comparison", "",
              "| Method | k=1 finite n | k=1 RMSE | k=2 finite n | k=2 RMSE |",
              "| --- | --- | --- | --- | --- |"]
    for row in result["conditional_rmse"]:
        a, b = row["first"], row["second"]
        def rmse(point):
            return "unavailable" if point["rmse"] is None else f"{point['rmse']:.4f}"
        lines.append(f"| {labels[row['method']]} | {a['available']} | {rmse(a)} | {b['available']} | {rmse(b)} |")
    lines += ["", "These RMSE values condition on method- and design-specific finite-estimate "
              "events. Missing boundary or optimizer estimates are retained as failures in "
              "availability, never assigned zero squared error. The values cannot establish "
              "an unconditional risk ranking. Squared-error tails are untrimmed; the "
              "[tail audit](tail_report.md) describes their influence. No cross-design RMSE "
              "significance claim or paired-record test is made.", "",
              "## Conservative Monte Carlo precision budget", "",
              "For the same eight-contrast family and equal sample sizes, the sufficient "
              "per-design count for an untruncated half width h is ceil(log(320)/h²). "
              "This distribution-free bound is conservative; these are proposed counts, "
              "not executed studies, minimum necessary counts or power calculations.", "",
              "| Target half width (pp) | Sufficient attempts per design |",
              "| --- | --- |"]
    lines += [f"| {100*row['half_width']:.0f} | {row['replications_per_design']:,} |"
              for row in result["precision_budget"]]
    lines += ["", "This analysis was specified after inspecting the earlier results and is "
              "not a newly preregistered confirmation. Its bounds cover this fixed eight-item "
              "family; they do not cover data-dependent selection among all past/future "
              "designs, repeated monitoring or all claims in the repository. The original "
              "inference algorithms, archived records and seeds are unchanged.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="comparison output folder")
    parser.add_argument("--check", action="store_true", help="reject absent/stale committed JSON or Markdown")
    args = parser.parse_args()
    try:
        result = build()
        outputs = {"design_comparison.json": json.dumps(result, indent=2, allow_nan=False) + "\n",
                   "design_comparison.md": render(result)}
        if args.check:
            for name, content in outputs.items():
                path = args.output / name
                if not path.exists():
                    raise ValueError(f"missing or stale report: {path}")
                if name.endswith(".json"):
                    try:
                        verifier.compare(verifier.read_json(path), result, "comparison")
                    except (ValueError, TypeError) as error:
                        raise ValueError(f"stale report: {path}: {error}") from error
                elif path.read_text(encoding="utf-8") != content:
                    raise ValueError(f"stale report: {path}")
        else:
            args.output.mkdir(parents=True, exist_ok=True)
            for name, content in outputs.items():
                (args.output / name).write_text(content, encoding="utf-8")
        print(f"PASS: {'checked' if args.check else 'wrote'} eight independent-design contrasts; "
              "2,000 archived records, no new simulations")
    except (ValueError, OSError, EOFError, KeyError, TypeError, OverflowError) as error:
        print(f"FAIL design comparison: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
