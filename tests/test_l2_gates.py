"""L2 gate sensitivity tests: normalisation order, sign convention and the rule.

Same principle as the L1 tests. Each case breaks one instruction from the
methods guide and requires the result to move. A test that only re-prints the
pipeline output would pass with a broken pipeline.

    python -m unittest discover -s tests -v
"""
from __future__ import annotations

import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), ".deps"))
os.environ.setdefault("MPLCONFIGDIR",
                      os.path.join(os.path.dirname(HERE), ".mplconfig"))

import l2  # noqa: E402

C = 0.03875489491553498
T = 300.0
DATA = os.path.join(os.path.dirname(HERE), "data", "raw", "regular_G18")

# A shortened gap grid keeps the suite fast; the decision rule does not depend
# on the step size.
FAST_GRID = np.array([0.50, 1.00, 1.25, 2.00, 2.50, 3.00, 4.25, 6.00, 8.00])


def library() -> dict:
    saved = l2.GAP_GRID_KT
    l2.GAP_GRID_KT = FAST_GRID
    try:
        return l2.build_library(C, T)
    finally:
        l2.GAP_GRID_KT = saved


def measured(joint: bool = False, absolute: bool = False) -> np.ndarray:
    raw = l2.load_repeats(os.path.join(DATA, "L2_data.csv"))
    rh = raw["RH"].mean(axis=1)
    rs = raw["RS"].mean(axis=1)
    if joint:
        scale = max(np.max(np.abs(rh)), np.max(rs))
        xy = np.column_stack([rh / scale, rs / scale])
    else:
        xy = l2.normalise_pair(rh, rs)
    if absolute:
        xy = np.column_stack([np.abs(xy[:, 0]), xy[:, 1]])
    return xy


class TestL2Gates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lib = library()
        cls.xy = measured()

    def test_zero_crossing_survives_normalisation(self):
        raw = l2.load_repeats(os.path.join(DATA, "L2_data.csv"))
        l2.gate_zero_crossing(raw)          # raises if the sign change is lost
        rh = self.xy[:, 0]
        self.assertTrue(np.any(rh > 0) and np.any(rh < 0))

    def test_joint_normalisation_changes_the_fit(self):
        """Normalising R_H and R_s together collapses the scale information."""
        good = l2.fit(self.xy, self.lib)
        bad = l2.fit(measured(joint=True), self.lib)
        self.assertGreater(abs(bad["r_gap"] - good["r_gap"]) / good["r_gap"], 0.10)

    def test_absolute_rh_loses_the_neutral_crossing(self):
        """|R_H| is not caught by the residual; it is caught by the trajectory.

        Folding the sign makes the loop fit marginally better (0.0148 against
        0.0159), so a residual-only check would wave the mistake through. What
        disappears is the assignment: E_F then never crosses zero, and the drift
        +73 -> -125 meV becomes +73 -> +128 meV.
        """
        good = l2.fit(self.xy, self.lib)
        bad = l2.fit(measured(absolute=True), self.lib)
        t_good = self.lib["EFG"][good["nearest"]] * 1e3
        t_bad = self.lib["EFG"][bad["nearest"]] * 1e3

        def crosses(t):
            return bool(np.any(np.sign(t[:-1]) * np.sign(t[1:]) < 0))

        self.assertTrue(crosses(t_good))
        self.assertFalse(crosses(t_bad))
        self.assertGreater(bad["r_gap"], 0.9 * good["r_gap"])

    def test_fifteen_percent_rule_rejects_a_gapless_curve(self):
        """Data lying on the cone must keep the cone even if a gap fits better."""
        cone = self.lib[("cone", 0.0)]
        pick = np.linspace(5, len(cone) - 6, 56).astype(int)
        rng = np.random.default_rng(7)
        on_cone = cone[pick] + rng.normal(0.0, 0.005, size=(56, 2))
        res = l2.fit(on_cone, self.lib)
        self.assertEqual(res["model"], "cone")
        self.assertGreaterEqual(res["r_gap"], 0.85 * res["r_cone"])

    def test_sign_convention_maps_positive_rh_to_positive_ef(self):
        """(n - p) is the course convention: R_H > 0 means E_F > 0."""
        pair = l2.model_pair(np.array([3.0, 1.0]), np.array([1.0, 3.0]))
        self.assertGreater(pair[0, 0], 0.0)
        self.assertLess(pair[1, 0], 0.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
