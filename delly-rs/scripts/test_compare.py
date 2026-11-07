#!/usr/bin/env python3
"""Regression tests for the comparison harness and C++ probe ordering.

Run with: python3 scripts/test_compare.py
The fake callers/reader isolate shell failure handling; real BCF comparisons are
exercised separately with scripts/compare.sh.
"""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class CompareTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="delly compare tests ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.out = self.root / "output files"
        self.out.mkdir()
        reader = self.root / "bcftools"
        reader.write_text("""#!/usr/bin/env python3
import pathlib, sys
data = pathlib.Path(sys.argv[-1]).read_bytes()
if data == b'invalid': sys.exit(1)
print('variant')
""")
        reader.chmod(0o755)
        self.env = dict(os.environ, PATH=str(self.root) + os.pathsep + os.environ["PATH"])

    def binary(self, name, mode="ok"):
        path = self.root / name
        path.write_text("""#!/usr/bin/env python3
from pathlib import Path
import sys
mode = MODE
if mode == 'fail': sys.exit(7)
args = sys.argv
out = Path(args[args.index('-o') + 1])
out.write_bytes(b'invalid' if mode == 'invalid' else b'calls')
if mode != 'missing-index':
    Path(str(out) + '.csi').write_bytes(b'different' if mode == 'index-diff' else b'index')
if '-d' in args:
    Path(args[args.index('-d') + 1]).write_bytes(b'different' if mode == 'dump-diff' else b'dump')
""".replace("MODE", repr(mode)))
        path.chmod(0o755)
        return path

    def run_compare(self, ref="ok", rs="ok", extra=()):
        return subprocess.run(["bash", str(ROOT / "scripts/compare.sh"),
            str(self.binary("reference caller", ref)), str(self.binary("rust caller", rs)),
            str(self.out), "test", *extra], env=self.env, capture_output=True, text=True)

    def test_failed_callers_cannot_pass_using_stale_outputs(self):
        for impl in ("ref", "rs"):
            for suffix in ("bcf", "bcf.csi"):
                (self.out / f"test.{impl}.{suffix}").write_bytes(b"stale")
        got = self.run_compare("fail", "fail")
        self.assertNotEqual(got.returncode, 0)
        self.assertNotIn("IDENTICAL", got.stdout)
        self.assertIn("exit 7", got.stderr)
        self.assertTrue((self.out / "test.rs.log").exists())

    def test_success_uses_fresh_files_and_supports_long_dump_option(self):
        (self.out / "test.ref.bcf").write_bytes(b"stale")
        got = self.run_compare(extra=("--dump=ignored.gz",))
        self.assertEqual(got.returncode, 0, got.stderr)
        self.assertIn("IDENTICAL", got.stdout)
        self.assertEqual((self.out / "test.ref.bcf").read_bytes(), b"calls")
        self.assertEqual((self.out / "test.rs.dump.gz").read_bytes(), b"dump")

    def test_missing_invalid_and_different_outputs_fail(self):
        for mode in ("fail", "missing-index", "invalid", "index-diff", "dump-diff"):
            with self.subTest(mode=mode):
                got = self.run_compare(rs=mode, extra=("-d", "ignored.gz"))
                self.assertNotEqual(got.returncode, 0, got.stdout)

    def test_output_override_is_rejected(self):
        got = self.run_compare(extra=("-o", str(self.root / "other.bcf")))
        self.assertEqual(got.returncode, 2)


class ProbeOrderTests(unittest.TestCase):
    def test_parallel_collection_retains_serial_sort_permutation(self):
        with tempfile.TemporaryDirectory(prefix="delly probe order ") as td:
            binary = Path(td) / "probe-order"
            subprocess.run(["c++", "-std=c++17", "-O2", str(ROOT / "tests/probe_order.cpp"),
                "-o", str(binary)], check=True)
            subprocess.run([str(binary)], check=True)


if __name__ == "__main__":
    unittest.main()
