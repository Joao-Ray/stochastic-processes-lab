"""Reject plausible archive corruption even when file checksums are refreshed.

These tests need only the standard library and the committed result bundle.
They exercise the end-to-end checker without importing the inference code.
"""
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("verify_noise_results", ROOT / "scripts/verify_noise_results.py")
VERIFIER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERIFIER)
BASE = ROOT / VERIFIER.STUDY


class TestResultIntegrity(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name) / "confirmation"
        shutil.copytree(BASE / "results/confirmation", self.folder)

    def verify(self):
        return VERIFIER.verify_bundle(self.folder, BASE / "confirmation.json", ROOT, require_extremes=True)

    def edit_json(self, name, change):
        path = self.folder / name
        data = json.loads(path.read_text())
        change(data)
        path.write_text(json.dumps(data) + "\n")

    def edit_records(self, change):
        path = self.folder / "replications.csv.gz"
        with gzip.open(path, "rt", newline="") as stream:
            reader = csv.DictReader(stream)
            names, rows = reader.fieldnames, list(reader)
        change(rows)
        with gzip.open(path, "wt", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=names)
            writer.writeheader()
            writer.writerows(rows)
        # A checksum alone cannot detect an internally inconsistent new archive.
        checksum = VERIFIER.digest(path)
        for file, field in (("summary.json", "record_archive_sha256"),
                            ("tails.json", "archive_sha256"), ("extreme_profiles.json", "archive_sha256")):
            self.edit_json(file, lambda data, key=field: data.update({key: checksum}))

    def test_both_committed_bundles_and_cli_receipt(self):
        main = VERIFIER.verify_bundle(BASE / "results/main", BASE / "config.json", ROOT)
        confirmation = self.verify()
        self.assertEqual((main["records"], main["planned_optimizer_audits"]), (7000, 56))
        self.assertEqual((confirmation["records"], confirmation["posthoc_extremes"]), (2000, 2))
        receipt = Path(self.temporary.name) / "receipt.json"
        process = subprocess.run([sys.executable, str(ROOT / "scripts/verify_noise_results.py"),
                                  "--output", str(receipt)], cwd=self.temporary.name,
                                 capture_output=True, text=True, timeout=30)
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(json.loads(receipt.read_text())["bundles"]["confirmation"], confirmation)

    def test_changed_rmse_rejected(self):
        self.edit_json("summary.json", lambda data: data["scenarios"][0]["methods"]["reml"]["point"].update(rmse=.1))
        with self.assertRaisesRegex(VERIFIER.IntegrityError, r"reml.point.rmse"):
            self.verify()

    def test_changed_wilson_interval_rejected(self):
        self.edit_json("summary.json", lambda data: data["scenarios"][0]["methods"]["reml"]["success_and_coverage"].update(wilson_95=[.99, 1.]))
        with self.assertRaisesRegex(VERIFIER.IntegrityError, "wilson_95"):
            self.verify()

    def test_changed_paired_population_rejected(self):
        self.edit_json("summary.json", lambda data: data["scenarios"][0]["paired"]["reml_minus_ml"].update(n=1000))
        with self.assertRaisesRegex(VERIFIER.IntegrityError, "paired.reml_minus_ml.n"):
            self.verify()

    def test_duplicate_replication_rejected_after_rehash(self):
        self.edit_records(lambda rows: rows[1].update(replication=rows[0]["replication"]))
        with self.assertRaisesRegex(VERIFIER.IntegrityError, "duplicate/missing replication"):
            self.verify()

    def test_missing_record_rejected_after_rehash(self):
        self.edit_records(lambda rows: rows.pop())
        with self.assertRaisesRegex(VERIFIER.IntegrityError, "record count"):
            self.verify()

    def test_unavailable_hit_rejected_after_rehash(self):
        self.edit_records(lambda rows: rows[0].update(reml_test_available="False"))
        with self.assertRaisesRegex(VERIFIER.IntegrityError, "hit from unavailable test"):
            self.verify()

    def test_hit_inconsistent_with_lr_rejected_after_rehash(self):
        self.edit_records(lambda rows: rows[0].update(reml_hit="False"))
        with self.assertRaisesRegex(VERIFIER.IntegrityError, "hit/cutoff"):
            self.verify()

    def test_nonfinite_estimate_rejected_after_rehash(self):
        self.edit_records(lambda rows: rows[0].update(reml_theta="NaN"))
        with self.assertRaisesRegex(VERIFIER.IntegrityError, "nonfinite field"):
            self.verify()

    def test_source_mismatch_rejected(self):
        self.edit_json("summary.json", lambda data: data["source_sha256"].update({"src/ou_measurement_error.py": "0" * 64}))
        with self.assertRaisesRegex(VERIFIER.IntegrityError, "source hash src/ou_measurement_error.py"):
            self.verify()

    def test_deleted_source_evidence_rejected(self):
        self.edit_json("summary.json", lambda data: data["source_sha256"].pop("src/ou_measurement_error.py"))
        with self.assertRaisesRegex(VERIFIER.IntegrityError, "required source list"):
            self.verify()

    def test_changed_tail_share_rejected(self):
        self.edit_json("tails.json", lambda data: data["scenarios"][0]["methods"]["reml"].update(top_five_squared_error_share=0.))
        with self.assertRaisesRegex(VERIFIER.IntegrityError, "top_five_squared_error_share"):
            self.verify()

    def test_deleted_planned_audit_rejected(self):
        self.edit_json("optimizer_audit.json", lambda data: data.pop())
        with self.assertRaisesRegex(VERIFIER.IntegrityError, "missing planned record"):
            self.verify()

    def test_false_optimizer_gap_rejected(self):
        self.edit_json("optimizer_audit.json", lambda data: data[0].update(gap=1.))
        with self.assertRaisesRegex(VERIFIER.IntegrityError, "audit.gap"):
            self.verify()

    def test_unresolved_extreme_gap_rejected(self):
        def change(data):
            data["examples"][0]["profile_lr"][0] = -.1
            data["examples"][0]["negative_profile_gap"] = True
        self.edit_json("extreme_profiles.json", change)
        with self.assertRaisesRegex(VERIFIER.IntegrityError, "unresolved numerical discrepancy"):
            self.verify()


if __name__ == "__main__":
    unittest.main()
