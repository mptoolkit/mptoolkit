#!/usr/bin/env python3
"""CLI regressions using test-only, deliberately noncanonical pheap fixtures."""

import argparse
import hashlib
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class CheckOrtho(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="mptk-check-ortho-")
        cls.root = Path(cls.temp.name)
        cls.env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
                       MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", MP_PRECISION="17")
        result = subprocess.run([str(ARGS.fixtures)], cwd=cls.root, env=cls.env,
                                capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def run_checker(self, name, *options, status=0):
        path = self.root / name
        before = hashlib.sha256(path.read_bytes()).digest()
        result = subprocess.run([str(ARGS.checker), *options, str(path)],
                                cwd=self.root, env=self.env, capture_output=True, text=True)
        self.assertEqual(result.returncode, status, result.stdout + result.stderr)
        self.assertEqual(hashlib.sha256(path.read_bytes()).digest(), before,
                         "checker modified its input")
        return result

    def metrics(self, name):
        result = self.run_checker(name, "--quiet")
        values = {key: (float(a), float(r)) for key, a, r in
                  (line.split() for line in result.stdout.splitlines())}
        self.assertEqual(set(values), {"max_left_orthogonality", "max_stored_density_bond",
                                       "right_density_fixed_point"})
        return values

    def close(self, actual, expected):
        self.assertTrue(math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12),
                        f"{actual} != {expected}")

    def test_canonical_and_symmetry_states(self):
        for name in ("canonical", "scaled", "u1-shift", "su2"):
            with self.subTest(name=name):
                for pair in self.metrics(name).values():
                    for value in pair:
                        self.close(value, 0.0)

    def test_tensor_error_is_not_repaired(self):
        for name in ("tensor-error", "su2-tensor-error"):
            with self.subTest(name=name):
                values = self.metrics(name)
                self.close(values["max_left_orthogonality"][0], 3 * math.sqrt(2))
                self.close(values["max_left_orthogonality"][1], 3)
                for key in ("max_stored_density_bond", "right_density_fixed_point"):
                    self.close(values[key][0], 3 / math.sqrt(2))
                    self.close(values[key][1], 3)

    def test_stale_interior_density_is_independent_of_orthogonality(self):
        values = self.metrics("stored-density-error")
        self.close(values["max_left_orthogonality"][0], 0)
        self.close(values["right_density_fixed_point"][0], 0)
        self.close(values["max_stored_density_bond"][0], 3 / math.sqrt(2))
        self.close(values["max_stored_density_bond"][1], 0.75)

    def test_wrong_boundary_density(self):
        values = self.metrics("fixed-point-error")
        self.close(values["max_left_orthogonality"][0], 0)
        for key in ("max_stored_density_bond", "right_density_fixed_point"):
            self.close(values[key][0], 0.4 * math.sqrt(2))
            self.close(values[key][1], 0.4 * math.sqrt(2) / math.sqrt(0.82))

    def test_zero_reference(self):
        values = self.metrics("zero-reference")
        self.close(values["max_stored_density_bond"][0], 1 / math.sqrt(2))
        self.assertEqual(values["max_stored_density_bond"][1], math.inf)

    def test_nonfinite_is_not_reported_as_zero(self):
        result = self.run_checker("nonfinite", "--quiet", status=1)
        self.assertIn("non-finite", result.stderr)

    def test_unsupported_representation(self):
        result = self.run_checker("right-canonical", status=1)
        self.assertIn("only InfiniteWavefunctionLeft", result.stderr)

    def test_verbose_output(self):
        result = self.run_checker("canonical", "-vv")
        for prefix in ("site 0 ", "site 1 ", "bond 0 ", "bond 1 ", "#quantity #abs #rel"):
            self.assertTrue(any(line.startswith(prefix) for line in result.stdout.splitlines()),
                            result.stdout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixtures", type=Path, required=True)
    parser.add_argument("--checker", type=Path, required=True)
    ARGS = parser.parse_args()
    ARGS.fixtures = ARGS.fixtures.resolve()
    ARGS.checker = ARGS.checker.resolve()
    unittest.main(argv=[sys.argv[0]], verbosity=2)
