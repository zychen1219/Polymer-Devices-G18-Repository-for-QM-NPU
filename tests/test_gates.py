"""Gate sensitivity tests.

A gate that never fires is a comment, not a check. Each test breaks one
physical assumption on purpose and requires the corresponding assertion to
raise. The tests do not touch the source files: the breakage is applied to the
module in memory and undone afterwards.

    python -m unittest discover -s tests -v
"""
from __future__ import annotations

import os
import sys
import unittest
from unittest import mock

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), ".deps"))
os.environ.setdefault("MPLCONFIGDIR",
                      os.path.join(os.path.dirname(HERE), ".mplconfig"))

import l1  # noqa: E402

PARAMS = {
    "c": 0.03875489491553498,
    "T": 300.0,
    "EF_op_eV": 0.17610827554768008,
    "rho_ohm_sq": 527.5077673801122,
}


def swapped_band(EF, c, T=300.0):
    """The wrong band-edge count: integrate from E_F upwards."""
    E = l1._E[l1._E >= EF]
    f = l1.fermi_dirac(E, EF, T)
    return l1.N_A * c * l1.trapz(E * f, E)


def swapped_fermi(EF, c, T=300.0):
    """The wrong Fermi count: integrate from the Dirac point upwards."""
    f = l1.fermi_dirac(l1._E, EF, T)
    return l1.N_A * c * l1.trapz(l1._E * f, l1._E)


class TestGates(unittest.TestCase):
    def test_reference_values_are_reproduced(self):
        for ef, target in [(0.0, 1.000), (0.073, 2.005),
                           (0.100, 2.599), (0.200, 5.104)]:
            self.assertAlmostEqual(l1.correction(ef, PARAMS["c"]), target, places=2)

    def test_swapped_limits_are_rejected(self):
        """Transposing the two integration limits makes C < 1 and must fail."""
        with mock.patch.object(l1, "n_band", swapped_band), \
                mock.patch.object(l1, "n_fermi", swapped_fermi):
            self.assertLess(l1.correction(0.073, PARAMS["c"]), 1.0)
            with self.assertRaises(AssertionError):
                l1.run(dict(PARAMS), [], verbose=False)

    def test_missing_dos_prefactor_is_rejected(self):
        """Dropping N_A shifts every count by 3.816e15 and must fail gate 3."""
        rows = [{"EF_eV": 0.05, "n_band_cm2": 1.0e12}]
        with mock.patch.object(l1, "N_A", 1.0):
            with self.assertRaises(AssertionError):
                l1.run(dict(PARAMS), rows, verbose=False)

    def test_working_point_consistency(self):
        res = l1.run(dict(PARAMS), [], verbose=False)
        self.assertGreater(res["correction_factor"], 1.0)
        self.assertAlmostEqual(res["mobility_corrected_cm2_Vs"]
                               / res["mobility_band_cm2_Vs"],
                               res["correction_factor"], places=6)


if __name__ == "__main__":
    unittest.main(verbosity=2)
