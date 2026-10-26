"""Independent-design contrasts must retain failures and Monte Carlo uncertainty."""
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("compare_noise_designs", ROOT / "scripts/compare_noise_designs.py")
ANALYSIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)


class TestDesignComparison(unittest.TestCase):
    def test_independent_contrast_unequal_samples(self):
        result = ANALYSIS.contrast(2, 10, 12, 20)
        self.assertAlmostEqual(result["difference"], .4)
        # Independent Bernoulli plug-in variance .16/10 + .24/20 = .028.
        self.assertAlmostEqual(result["mc_se"], .1673320053068151)
        self.assertEqual((result["first"]["n"], result["second"]["n"]), (10, 20))

    def test_boundary_counts_keep_nonzero_interval(self):
        result = ANALYSIS.contrast(1000, 1000, 1000, 1000)
        self.assertEqual(result["mc_se"], 0.)
        # A zero plug-in SE cannot make the uncertainty interval collapse.
        self.assertAlmostEqual(result["simultaneous_interval"][0], -.07594946343321836)
        self.assertAlmostEqual(result["simultaneous_interval"][1], .07594946343321836)
        clipped = ANALYSIS.contrast(0, 1, 1, 1)["simultaneous_interval"]
        self.assertEqual(clipped, [-1., 1.])

    def test_reversing_designs_reverses_interval(self):
        forward = ANALYSIS.contrast(7, 20, 15, 30)
        backward = ANALYSIS.contrast(15, 30, 7, 20)
        self.assertAlmostEqual(forward["difference"], -backward["difference"])
        lo, hi = forward["simultaneous_interval"]
        self.assertEqual(backward["simultaneous_interval"], [-hi, -lo])

    def test_simultaneous_family_widens_bound(self):
        one = ANALYSIS.contrast(20, 40, 20, 40, family=1)
        eight = ANALYSIS.contrast(20, 40, 20, 40, family=8)
        self.assertGreater(eight["half_width"], one["half_width"])

    def test_invalid_inputs_rejected(self):
        for counts in ((-1, 10, 2, 10), (11, 10, 2, 10), (0, 0, 1, 10),
                       (1., 10, 1, 10), (True, 10, 1, 10)):
            with self.subTest(counts=counts), self.assertRaises(ValueError):
                ANALYSIS.contrast(*counts)
        for options in ({"alpha": 0}, {"alpha": 1}, {"alpha": float("nan")},
                        {"family": 0}, {"family": 1.5}, {"family": True}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                ANALYSIS.contrast(1, 10, 2, 10, **options)

    def test_small_binomial_coverage_by_exact_enumeration(self):
        # Enumerate independent Bin(12,p1), Bin(9,p2), no Monte Carlo randomness.
        # This catches an under-scaled radius even with skewed outcomes.
        for p1, p2 in ((.1, .1), (.1, .8), (.5, .5), (.95, .2)):
            failure = 0.
            for x in range(13):
                for y in range(10):
                    interval = ANALYSIS.contrast(x, 12, y, 9, alpha=.2, family=1)["simultaneous_interval"]
                    if not interval[0] <= p2 - p1 <= interval[1]:
                        failure += (math.comb(12, x) * p1**x * (1-p1)**(12-x)
                                    * math.comb(9, y) * p2**y * (1-p2)**(9-y))
            self.assertLessEqual(failure, .2 + 1e-12)

    def test_precision_budget_rounds_up(self):
        self.assertEqual(ANALYSIS.required_per_design(.02), 14421)
        for value in (0, -.1, 1.1, float("nan")):
            with self.assertRaises(ValueError):
                ANALYSIS.required_per_design(value)

    def test_archived_confirmation_counts_and_conditional_risk(self):
        result = ANALYSIS.build()
        self.assertEqual(len(result["contrasts"]), 8)
        coverage = next(x for x in result["contrasts"] if x["method"] == "reml"
                        and x["endpoint"] == "success_and_coverage")
        self.assertEqual((coverage["first"]["hits"], coverage["second"]["hits"]), (957, 946))
        self.assertAlmostEqual(coverage["difference"], -.011)
        availability = next(x for x in result["contrasts"] if x["method"] == "reml"
                            and x["endpoint"] == "finite_point_available")
        self.assertEqual((availability["first"]["hits"], availability["second"]["hits"]), (973, 983))
        reml = next(x for x in result["conditional_rmse"] if x["method"] == "reml")
        self.assertEqual((reml["first"]["available"], reml["second"]["available"]), (973, 983))
        self.assertAlmostEqual(reml["first"]["rmse"], 1.5394, places=4)
        self.assertAlmostEqual(reml["second"]["rmse"], 1.2986, places=4)

    def test_cli_generation_and_stale_report_rejection(self):
        with tempfile.TemporaryDirectory() as temporary:
            command = [sys.executable, str(ROOT / "scripts/compare_noise_designs.py"), "--output", temporary]
            process = subprocess.run(command, cwd=temporary, capture_output=True, text=True, timeout=30)
            self.assertEqual(process.returncode, 0, process.stderr)
            process = subprocess.run(command + ["--check"], capture_output=True, text=True, timeout=30)
            self.assertEqual(process.returncode, 0, process.stderr)
            path = Path(temporary) / "design_comparison.md"
            path.write_text(path.read_text().replace("1.5394", "0.0000"))
            process = subprocess.run(command + ["--check"], capture_output=True, text=True, timeout=30)
            self.assertEqual(process.returncode, 1)
            self.assertIn("stale", process.stderr)

    def test_json_check_allows_roundoff_but_rejects_changed_count(self):
        with tempfile.TemporaryDirectory() as temporary:
            command = [sys.executable, str(ROOT / "scripts/compare_noise_designs.py"), "--output", temporary]
            subprocess.run(command, capture_output=True, check=True, timeout=30)
            path = Path(temporary) / "design_comparison.json"
            data = json.loads(path.read_text())
            # Different libm versions may differ in the final floating-point bit.
            data["contrasts"][0]["half_width"] += 1e-15
            path.write_text(json.dumps(data))
            process = subprocess.run(command + ["--check"], capture_output=True, text=True, timeout=30)
            self.assertEqual(process.returncode, 0, process.stderr)
            data["contrasts"][0]["first"]["hits"] += 1
            path.write_text(json.dumps(data))
            process = subprocess.run(command + ["--check"], capture_output=True, text=True, timeout=30)
            self.assertEqual(process.returncode, 1)
            self.assertIn("stale", process.stderr)


if __name__ == "__main__":
    unittest.main()
