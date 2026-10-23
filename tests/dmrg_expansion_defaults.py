#!/usr/bin/env python3
"""Check the effective DMRG expansion defaults and explicit overrides."""

import argparse
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest


class ExpansionDefaults(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="mptk-dmrg-expansion-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
                        MKL_NUM_THREADS="1", MP_PRECISION="17")
        self.run_tool(ARGS.model, "-o", "lat")
        self.serial = 0

    def run_tool(self, program, *args):
        result = subprocess.run([str(program), *map(str, args)], cwd=self.root,
                                env=self.env, text=True, capture_output=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result.stdout

    def dmrg(self, *options):
        self.serial += 1
        state = f"psi-{self.serial}"
        self.run_tool(ARGS.construct, "-l", "lat", "-o", state, "--finite", "0.5:-0.5:0.5:-0.5")
        return self.run_tool(ARGS.dmrg, "-w", state, "-H", "lat:H_J1", "-m", "4x2", *options)

    def check_settings(self, output, side, algorithm, factor, sector):
        match = re.search(rf"^Using {side}-expansion algorithm: {algorithm} with "
                          r"increment: (\S+) expansion factor: (\S+) per sector: (\d+)$",
                          output, re.MULTILINE)
        self.assertIsNotNone(match, output)
        self.assertEqual(float(match[1]), 1.0)
        self.assertAlmostEqual(float(match[2]), float(factor), delta=1e-14)
        self.assertEqual(int(match[3]), sector)

    def test_default_post_factor(self):
        self.check_settings(self.dmrg(), "post", "rsvd", "0.1", 1)

    def test_algorithm_sector_defaults(self):
        for algorithm, sector in (("fullsvd", 1), ("rsvd", 1), ("range", 0), ("random", 0)):
            with self.subTest(algorithm=algorithm):
                output = self.dmrg("--pre", algorithm, "--post", algorithm)
                for side in ("pre", "post"):
                    self.check_settings(output, side, algorithm, "0.1", sector)

    def test_explicit_settings_override_defaults(self):
        for algorithm, sector in (("fullsvd", 0), ("rsvd", 0), ("range", 1), ("random", 1)):
            with self.subTest(algorithm=algorithm):
                output = self.dmrg("--pre", algorithm, "--post", algorithm,
                                   "--pre-sector", sector, "--post-sector", sector,
                                   "--pre-factor", "0.25", "--post-factor", "0.25")
                for side in ("pre", "post"):
                    self.check_settings(output, side, algorithm, "0.25", sector)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for tool in ("model", "construct", "dmrg"):
        parser.add_argument("--" + tool, type=Path, required=True)
    ARGS, remaining = parser.parse_known_args()
    ARGS.model = ARGS.model.resolve()
    ARGS.construct = ARGS.construct.resolve()
    ARGS.dmrg = ARGS.dmrg.resolve()
    unittest.main(argv=[__file__, *remaining])
