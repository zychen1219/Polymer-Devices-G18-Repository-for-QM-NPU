"""Assertions for the L3 inversion.

Each test pins one physical statement that the reconstruction rests on, so a
regression shows up as a failed assertion rather than a plausible-looking map.

    python -m unittest discover -s tests -v
"""
import json
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import l3  # noqa: E402

DATA = os.path.join(ROOT, "data", "raw", "regular_G18")


def load():
    A = np.load(os.path.join(DATA, "L3_A.npy"))
    D = np.loadtxt(os.path.join(DATA, "L3_d.csv"), delimiter=",")
    meta = json.load(open(os.path.join(DATA, "L3_meta.json"), encoding="utf-8"))
    return A, D, meta


class TestForwardMatrix(unittest.TestCase):
    def test_rows_sum_to_zero(self):
        """A m is blind to a constant in m, so the map's mean is unidentifiable."""
        A, _, _ = load()
        self.assertLess(np.abs(A.sum(axis=1)).max(), 1e-11)

    def test_data_mean_is_not_removed(self):
        """The zero-mean rule applies to m, not to d.

        A's rows sum to zero, but its columns do not, so A m carries the level of
        the data.  Centring d before the fit forces the model to chase a target no
        map can produce and raises the residual by two orders of magnitude.
        """
        A, D, _ = load()
        L = l3.gradient_matrix(10)
        dbar = D.mean(axis=0)
        eta = l3.eta0(A, L)
        m_raw = l3.solve(A, L, dbar, eta)
        m_cen = l3.solve(A, L, dbar - dbar.mean(), eta)
        r_raw = np.linalg.norm(A @ m_raw - dbar)
        r_cen = np.linalg.norm(A @ m_cen - (dbar - dbar.mean()))
        self.assertLess(r_raw, r_cen / 10.0)


class TestRegularisation(unittest.TestCase):
    def test_solution_has_zero_mean(self):
        A, D, _ = load()
        m = l3.solve(A, l3.gradient_matrix(10), D.mean(axis=0), 1e-3)
        self.assertLess(abs(m.mean()), 1e-12)

    def test_maximum_shrinks_with_eta(self):
        """More smoothing can only lower the peak of the reconstruction."""
        A, D, _ = load()
        L = l3.gradient_matrix(10)
        dbar = D.mean(axis=0)
        peaks = [np.abs(l3.solve(A, L, dbar, f * l3.eta0(A, L))).max()
                 for f in (1.0, 10.0, 100.0)]
        self.assertTrue(peaks[0] >= peaks[1] >= peaks[2])
        self.assertLess(peaks[2], 0.9 * peaks[0])

    def test_decision_is_stable_to_a_decade_of_eta(self):
        """The verdict must not depend on a hand-tuned parameter."""
        A, D, _ = load()
        L = l3.gradient_matrix(10)
        dbar = D.mean(axis=0)
        base = l3.eta0(A, L)
        for f in (0.1, 1.0, 10.0):
            m = l3.solve(A, L, dbar, f * base)
            self.assertGreater(np.abs(m).max(), l3.TAU,
                               f"verdict flips at {f} eta0")


class TestResolution(unittest.TestCase):
    def test_gaussian_profile_recovers_its_width(self):
        """A profile built with FWHM 2.0 px must be measured as such."""
        x = np.arange(11.0)
        profile = np.exp(-4.0 * np.log(2.0) * (x - 5.0) ** 2 / 2.0 ** 2)
        got = l3.fwhm_profile(profile, peak=5)
        self.assertAlmostEqual(got["interp"], 2.0, delta=0.05)
        self.assertAlmostEqual(got["gauss"], 2.0, delta=0.10)

    def test_integer_count_underestimates(self):
        """The peak sits between samples, so counting samples gives 1 px."""
        A, D, _ = load()
        L = l3.gradient_matrix(10)
        eta = l3.eta0(A, L)
        res = l3.resolution(A, L, eta)
        for item in res.values():
            self.assertLess(item["coarse"], item["interp"])
            self.assertGreater(item["interp"], 1.0)


class TestDetectionFloor(unittest.TestCase):
    def test_reproducible_and_positive(self):
        A, D, _ = load()
        L = l3.gradient_matrix(10)
        eta = l3.eta0(A, L)
        a = l3.detection_limit(A, L, eta, D, n=20, seed=7)
        b = l3.detection_limit(A, L, eta, D, n=20, seed=7)
        self.assertEqual(a["delta95"], b["delta95"])
        self.assertGreater(a["delta95"], 0.0)

    def test_noise_floor_is_below_signal(self):
        """The map's extreme must sit above the floor it is judged against."""
        A, D, _ = load()
        L = l3.gradient_matrix(10)
        eta = l3.eta0(A, L)
        dbar = D.mean(axis=0)
        det = l3.detection_limit(A, L, eta, D, n=20, seed=7)
        m = l3.solve(A, L, dbar, eta)
        self.assertGreater(np.abs(m).max() / det["delta95"], 10.0)


if __name__ == "__main__":
    unittest.main()
