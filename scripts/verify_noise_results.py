"""Independently check archived noisy-OU results using only the Python stdlib.

This verifies records and reported summaries, not the statistical model or the
global optimum of a likelihood. No inference module or report generator is used.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
STUDY = Path("subprojects/ou-inference-study/measurement-error")
METHODS = ("conditional", "stationary_ignorant", "ml", "reml")
PAIRS = (("ml", "conditional"), ("reml", "conditional"), ("reml", "ml"))
SOURCES = {
    (STUDY / name).as_posix() for name in ("run.py", "THEORY.md", "PROTOCOL.md", "AMENDMENTS.md")
} | {f"src/{name}.py" for name in
     ("ou_measurement_error", "ou_inference", "ou_uncertainty", "ou_monte_carlo")}


class IntegrityError(ValueError):
    """A check failed, with a location suitable for a CI diagnostic."""


def require(condition, where):
    if not condition:
        raise IntegrityError(where)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    def invalid(value):
        raise IntegrityError(f"{path}: nonfinite JSON value {value}")
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=invalid)


def compare(actual, expected, where):
    """Counts/types are exact; floating summaries allow summation-order error."""
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and actual.keys() == expected.keys(), where + ": keys")
        for key, value in expected.items():
            compare(actual[key], value, f"{where}.{key}")
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(actual) == len(expected), where + ": length")
        for i, value in enumerate(expected):
            compare(actual[i], value, f"{where}[{i}]")
    elif isinstance(expected, bool) or expected is None or isinstance(expected, str):
        require(type(actual) is type(expected) and actual == expected, where)
    elif isinstance(expected, int):
        require(type(actual) is int and actual == expected, where)
    else:
        require(type(actual) in (int, float) and math.isfinite(actual)
                and math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-11), where)


def number(row, field):
    text = row[field]
    if text == "":
        return None
    try:
        value = float(text)
    except (ValueError, TypeError) as error:
        raise IntegrityError(f"invalid numeric field {field}: {text!r}") from error
    require(math.isfinite(value), f"nonfinite field {field}")
    return value


def boolean(row, field):
    require(row[field] in ("True", "False"), f"invalid boolean field {field}")
    return row[field] == "True"


def proportion(hits):
    n = len(hits)
    if not n:
        return dict(n=0, value=None, mc_se=None, wilson_95=None)
    p, z = sum(hits) / n, 1.959963984540054
    denominator = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denominator
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denominator
    return dict(n=n, value=p, mc_se=math.sqrt(p * (1 - p) / n),
                wilson_95=[center - half, center + half])


def point_metrics(values, truth):
    errors = [v - truth for v in values if v is not None]
    n = len(errors)
    if not n:
        return dict(available=0, bias=None, rmse=None, bias_mc_se=None, rmse_mc_se=None)
    squared = [e * e for e in errors]
    rmse = math.sqrt(statistics.fmean(squared))
    return dict(available=n, bias=statistics.fmean(errors), rmse=rmse,
                bias_mc_se=statistics.stdev(errors) / math.sqrt(n) if n > 1 else None,
                rmse_mc_se=statistics.stdev(squared) / (2 * math.sqrt(n) * rmse)
                if n > 1 and rmse > 0 else None)


def source_hashes(recorded, expected, root, where):
    require(set(recorded) == expected, where + ": required source list")
    for relative, expected_digest in recorded.items():
        path = (root / relative).resolve()
        require(path.is_relative_to(root.resolve()), where + ": source outside repository")
        require(digest(path) == expected_digest, f"{where}: source hash {relative}")


def check_audit(audit, row, method, seed, where):
    compare(audit["seed"], seed, where + ".seed")
    compare(audit["objective_local"], number(row, method + "_nll"), where + ".objective_local")
    alternative = audit["objective_alternative"]
    require(type(alternative) in (int, float) and math.isfinite(alternative), where + ": alternative")
    gap = max(0., audit["objective_local"] - alternative)
    compare(audit["gap"], gap, where + ".gap")
    compare(audit["discrepancy"], gap > 1e-4, where + ".discrepancy")
    require(type(audit["alternative_converged"]) is bool, where + ": convergence type")
    require(0 <= audit["alternative_phi"] < 1 and 0 <= audit["alternative_gamma"] <= 1,
            where + ": alternative parameters")


def check_row(row, dt, repeats, level):
    # Chi-square(1) is the square of a standard normal: an independent cutoff.
    cutoff = statistics.NormalDist().inv_cdf((1 + level) / 2) ** 2
    for method in METHODS:
        point = number(row, method + "_theta")
        require(point is None or point > 0, method + ": invalid drift")
        available = boolean(row, method + "_test_available")
        hit = boolean(row, method + "_hit")
        require(not hit or available, method + ": hit from unavailable test")
        if method not in ("ml", "reml"):
            continue
        flags = {key: boolean(row, method + "_" + key) for key in
                 ("white_limit", "search_cap", "noise_boundary", "optimizer_converged", "optimizer_gap")}
        fields = {key: number(row, method + "_" + key) for key in
                  ("phi", "gamma", "nll", "mu", "scale", "lr", "noise_variance", "latent_variance")}
        exception = bool(row[method + "_exception"])
        if exception:
            require(all(value is None for value in fields.values()) and point is None
                    and not available and not flags["optimizer_converged"], method + ": exception fields")
            continue
        require(all(fields[key] is not None for key in ("phi", "gamma", "nll", "mu", "scale")),
                method + ": missing fit")
        require(0 <= fields["phi"] < 1 and 0 <= fields["gamma"] <= 1 and fields["scale"] > 0,
                method + ": invalid shape or scale")
        compare(flags["white_limit"], fields["phi"] == 0 or fields["gamma"] <= 1e-8,
                method + ": white limit")
        identified = flags["optimizer_converged"] and not flags["white_limit"] and not flags["search_cap"]
        require((point is not None) == identified, method + ": drift availability")
        if point is not None:
            compare(point, -math.log(fields["phi"]) / dt, method + ": drift/phi")
        if flags["white_limit"] and repeats == 1:
            require(fields["latent_variance"] is None and fields["noise_variance"] is None,
                    method + ": unidentified variance decomposition")
        else:
            compare(fields["latent_variance"], fields["gamma"] * fields["scale"], method + ": latent variance")
            compare(fields["noise_variance"], (1 - fields["gamma"]) * fields["scale"], method + ": noise variance")
        require(not available or (flags["optimizer_converged"] and not flags["optimizer_gap"]),
                method + ": membership/optimizer status")
        require((fields["lr"] is not None) == available, method + ": LR availability")
        require(fields["lr"] is None or fields["lr"] >= 0, method + ": negative available LR")
        compare(hit, available and fields["lr"] <= cutoff, method + ": hit/cutoff")


def check_metrics(block, scenario, truth, where):
    for method in METHODS:
        prefix = method + "_"
        actual = scenario["methods"][method]
        values = [number(row, prefix + "theta") for row in block]
        hits = [boolean(row, prefix + "hit") for row in block]
        available_hits = [boolean(row, prefix + "hit") for row in block
                          if boolean(row, prefix + "test_available")]
        compare(actual["point"], point_metrics(values, truth), where + "." + method + ".point")
        compare(actual["test_available"], len(available_hits), where + "." + method + ".test_available")
        compare(actual["conditional_coverage"], proportion(available_hits), where + "." + method + ".conditional_coverage")
        compare(actual["success_and_coverage"], proportion(hits), where + "." + method + ".all_attempts")
        if method in ("ml", "reml"):
            diagnostics = {flag: sum(boolean(row, prefix + flag) for row in block)
                           for flag in ("white_limit", "search_cap", "noise_boundary", "optimizer_gap")}
            diagnostics["exception"] = sum(bool(row[prefix + "exception"]) for row in block)
            diagnostics["nonconvergence"] = sum(not boolean(row, prefix + "optimizer_converged") for row in block)
            compare(actual["diagnostics"], diagnostics, where + "." + method + ".diagnostics")
            noise = [number(row, prefix + "noise_variance") for row in block]
            noise = [v for v in noise if v is not None]
            compare(actual["noise_variance_n"], len(noise), where + "." + method + ".noise_variance_n")
            compare(actual["noise_variance_mean"], statistics.fmean(noise) if noise else None,
                    where + "." + method + ".noise_variance_mean")
    require(set(scenario["methods"]) == set(METHODS), where + ": methods")
    require(set(scenario["paired"]) == {a + "_minus_" + b for a, b in PAIRS}, where + ": pairs")
    for a, b in PAIRS:
        paired = [(number(row, a + "_theta"), number(row, b + "_theta")) for row in block]
        differences = [(x - truth) ** 2 - (y - truth) ** 2 for x, y in paired
                       if x is not None and y is not None]
        n = len(differences)
        compare(scenario["paired"][a + "_minus_" + b],
                dict(n=n, mse_difference=statistics.fmean(differences) if n else None,
                     mc_se=statistics.stdev(differences) / math.sqrt(n) if n > 1 else None),
                where + ".paired." + a + "_minus_" + b)


def quantile(sorted_values, p):
    position = (len(sorted_values) - 1) * p
    left = math.floor(position)
    fraction = position - left
    return sorted_values[left] * (1 - fraction) + sorted_values[math.ceil(position)] * fraction


def check_tails(folder, blocks, truth, archive_sha, root):
    tails = read_json(folder / "tails.json")
    compare(tails["archive_sha256"], archive_sha, "tails.archive_sha256")
    compare(tails["source_sha256"], digest(root / STUDY / "tail_analysis.py"), "tails.source_sha256")
    require(len(tails["scenarios"]) == len(blocks), "tails: scenarios")
    seen = set()
    for scenario in tails["scenarios"]:
        name = scenario["id"]
        require(name in blocks and name not in seen, "tails: duplicate/unknown scenario")
        seen.add(name)
        require(set(scenario["methods"]) == set(METHODS), "tails: methods")
        for method in METHODS:
            values = [(int(row["replication"]), number(row, method + "_theta")) for row in blocks[name]]
            values = [(rep, v, abs(v - truth)) for rep, v in values if v is not None]
            where = "tails." + name + "." + method
            actual = scenario["methods"][method]
            if not values:
                compare(actual, dict(available=0), where)
                continue
            errors = sorted(error for _, _, error in values)
            squared_sum = math.fsum(error ** 2 for error in errors)
            largest = sorted(values, key=lambda value: value[2], reverse=True)[:5]
            expected = dict(available=len(values), median_absolute_error=statistics.median(errors),
                            q99_absolute_error=quantile(errors, .99), max_absolute_error=errors[-1],
                            top_five_squared_error_share=math.fsum(e ** 2 for _, _, e in largest) / squared_sum
                            if squared_sum else 0.)
            for key, value in expected.items():
                compare(actual[key], value, where + "." + key)
            # Ties may use a different argsort order; check membership, errors and ranks.
            selected = actual["largest_records"]
            require(len(selected) == len(largest), where + ": largest record count")
            lookup = {rep: (v, error) for rep, v, error in values}
            require(len({item["replication"] for item in selected}) == len(selected), where + ": duplicate tail record")
            for item, (_, _, ranked_error) in zip(selected, largest):
                rep = item["replication"]
                require(rep in lookup, where + ": unknown tail record")
                compare(item["theta"], lookup[rep][0], where + ".theta")
                compare(item["absolute_error"], lookup[rep][1], where + ".absolute_error")
                compare(item["absolute_error"], ranked_error, where + ".rank")


def check_extremes(folder, blocks, cfg, archive_sha, root):
    result = read_json(folder / "extreme_profiles.json")
    compare(result["posthoc"], True, "extremes.posthoc")
    compare(result["archive_sha256"], archive_sha, "extremes.archive_sha256")
    compare(result["config_sha256"], digest(root / STUDY / "confirmation.json"), "extremes.config_sha256")
    source_hashes(result["source_sha256"], {(STUDY / name).as_posix() for name in ("run.py", "diagnose_extremes.py")}
                  | {"src/ou_measurement_error.py"}, root, "extremes")
    require(len(result["examples"]) == len(blocks), "extremes: scenario count")
    seen = set()
    for example in result["examples"]:
        name = example["scenario"]
        require(name in blocks and name not in seen, "extremes: duplicate/unknown scenario")
        seen.add(name)
        index = list(blocks).index(name)
        design = cfg["scenarios"][index]
        selected = [row for row in blocks[name] if row["reml_theta"] != ""]
        worst = max(selected, key=lambda row: abs(number(row, "reml_theta") - cfg["theta"]))
        compare(example["replication"], int(worst["replication"]), "extremes.replication")
        compare(example["theta_estimate"], number(worst, "reml_theta"), "extremes.theta")
        compare(example["truth_theta"], float(cfg["theta"]), "extremes.truth")
        sites = cfg["budget"] // design["repeats"]
        compare(example["dt"], design["duration"] / (sites - 1), "extremes.dt")
        require(len(example["raw_observations"]) == sites, "extremes: sites")
        for values in example["raw_observations"]:
            require(len(values) == design["repeats"] and all(math.isfinite(x) for x in values), "extremes: data")
        check_audit(example["audit"], worst, "reml", cfg["seed"] + 950000 + index, "extremes.audit")
        grid, lr = example["theta_grid"], example["profile_lr"]
        require(len(grid) == len(lr) and len(grid) >= 151, "extremes: grid size")
        require(all(math.isfinite(x) and x > 0 for x in grid)
                and all(x < y for x, y in zip(grid, grid[1:]))
                and all(math.isfinite(x) for x in lr), "extremes: invalid profile grid")
        compare(example["negative_profile_gap"], min(lr) < -1e-6, "extremes.negative_profile_gap")
        require(not example["negative_profile_gap"] and not example["audit"]["discrepancy"],
                "extremes: unresolved numerical discrepancy")
    return len(result["examples"])


def verify_bundle(folder, config_path, root=ROOT, require_extremes=False):
    folder, root = Path(folder), Path(root)
    cfg, summary = read_json(Path(config_path)), read_json(folder / "summary.json")
    compare(summary["config"], cfg, "summary.config")
    compare(summary["config_sha256"], digest(Path(config_path)), "config hash")
    source_hashes(summary["source_sha256"], SOURCES, root, "summary")
    archive = folder / "replications.csv.gz"
    archive_sha = digest(archive)
    compare(summary["record_archive_sha256"], archive_sha, "archive hash")
    with gzip.open(archive, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames is not None and len(set(reader.fieldnames)) == len(reader.fieldnames),
                "archive: duplicate or missing columns")
        rows = list(reader)
    require(all(None not in row and None not in row.values() for row in rows), "archive: malformed CSV")
    compare(len(rows), len(cfg["scenarios"]) * cfg["replications"], "archive: record count")
    require(len(summary["scenarios"]) == len(cfg["scenarios"]), "summary: scenario count")
    blocks, by_record, expected_audits = {}, {}, set()
    for index, (design, scenario) in enumerate(zip(cfg["scenarios"], summary["scenarios"])):
        name = f"T={design['duration']:g},k={design['repeats']},eta={design['noise_sd']:g},stress={design['stress']}"
        compare(scenario["id"], name, "scenario id")
        require(name not in blocks, "duplicate design")
        for key, value in design.items():
            compare(scenario[key], value, name + "." + key)
        sites = cfg["budget"] // design["repeats"]
        dt = design["duration"] / (sites - 1)
        compare(scenario["sites"], sites, name + ".sites")
        compare(scenario["dt"], dt, name + ".dt")
        compare(scenario["attempted"], cfg["replications"], name + ".attempted")
        compare(scenario["model_assumptions_met"], design["stress"] == "none", name + ".assumptions")
        block = [row for row in rows if row["scenario"] == name]
        require(len(block) == cfg["replications"], name + ": missing records")
        ids = [int(row["replication"]) for row in block]
        require(sorted(ids) == list(range(cfg["replications"])), name + ": duplicate/missing replication")
        blocks[name] = sorted(block, key=lambda row: int(row["replication"]))
        for row in block:
            check_row(row, dt, design["repeats"], cfg["level"])
            rep = int(row["replication"])
            by_record[name, rep] = row
            if rep < cfg["audit_records_per_scenario"]:
                expected_audits.update((name, rep, method) for method in ("ml", "reml")
                                       if row[method + "_nll"] != "")
        check_metrics(block, scenario, cfg["theta"], name)
    require(set(row["scenario"] for row in rows) == set(blocks), "archive: unknown scenario")
    audits = read_json(folder / "optimizer_audit.json")
    seen = set()
    for audit in audits:
        key = audit["scenario"], audit["replication"], audit["method"]
        require(key in expected_audits and key not in seen, "audit: duplicate/unexpected record")
        seen.add(key)
        name, rep, method = key
        seed = cfg["seed"] + 10000 + list(blocks).index(name) * 20 + rep * 2 + int(method == "reml")
        check_audit(audit, by_record[name, rep], method, seed, "audit")
    require(seen == expected_audits, "audit: missing planned record")
    discrepancies = sum(audit["discrepancy"] for audit in audits)
    compare(summary["optimizer_audit_discrepancies"], discrepancies, "summary.audit_discrepancies")
    require(discrepancies == 0, "audit: unresolved numerical discrepancy")
    check_tails(folder, blocks, cfg["theta"], archive_sha, root)
    extremes = check_extremes(folder, blocks, cfg, archive_sha, root) if require_extremes else 0
    return dict(records=len(rows), scenarios=len(blocks), planned_optimizer_audits=len(audits),
                posthoc_extremes=extremes, archive_sha256=archive_sha,
                config_sha256=summary["config_sha256"], status="pass")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT, help="repository root")
    parser.add_argument("--output", type=Path, help="optional JSON validation receipt")
    args = parser.parse_args()
    try:
        base = args.root / STUDY
        result = dict(validator_sha256=digest(Path(__file__)), bundles={})
        for kind, config in (("main", "config.json"), ("confirmation", "confirmation.json")):
            result["bundles"][kind] = verify_bundle(base / "results" / kind, base / config,
                                                  args.root, require_extremes=kind == "confirmation")
            checked = result["bundles"][kind]
            print(f"PASS {kind}: {checked['records']} records, {checked['scenarios']} designs, "
                  f"{checked['planned_optimizer_audits']} planned audits, {checked['posthoc_extremes']} post hoc examples")
        print("Checks cover archive consistency and provenance; they do not certify model assumptions, "
              "global likelihood optima, originality or publication readiness.")
        if args.output:
            args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    except (IntegrityError, OSError, EOFError, KeyError, TypeError, ValueError, OverflowError) as error:
        print(f"FAIL archive integrity: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
