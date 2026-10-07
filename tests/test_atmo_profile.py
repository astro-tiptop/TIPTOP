import unittest
import warnings
import numpy as np

from tiptop.atmoProfile import (
    GRID, H_GL, F_GL_MIN, F_GL_MEDIAN, F_GL_RANGE, DEFAULT_LAYERS,
    logGrid, layerMatrix, seeingToR0, r0ToSeeing, theta0FromProfile,
    generateProfile, layerWind,
)

"""
Unit tests for tiptop.atmoProfile (parametric normalized Cn2 profiles).
theta0 and the GL fraction are imposed exactly by construction, so they are
checked with tight tolerances on several layer sets and random draws.
"""

# 5 conditions of Agapito+ 2021: (seeing [arcsec], theta0 [arcsec], GL fraction)
CONDITIONS = [
    (0.509, 2.82, 0.766),
    (0.603, 2.43, 0.690),
    (0.725, 2.00, 0.698),
    (0.889, 1.62, 0.534),
    (1.05, 1.35, 0.291),
]
LAYER_SETS = {
    'default': DEFAULT_LAYERS,
    'grid': None,
    'coarse': [0, 400, 2000, 5000, 9000, 13000, 17000, 21000],
}
RTOL = 1e-10


def _glFraction(heights, weights) -> float:
    """Sum of the weights at heights <= H_GL."""
    return float(np.sum(np.asarray(weights)[np.asarray(heights) <= H_GL]))


def _rngs():
    """None (median profile) and a few seeded generators."""
    return [None] + [np.random.default_rng(s) for s in (1, 2, 3)]


class TestGenerateProfileExactness(unittest.TestCase):

    def _check(self, useGl, glOnlyChecks=True):
        for seeing, theta0, gl in CONDITIONS:
            r0 = seeingToR0(seeing)
            for name, layers in LAYER_SETS.items():
                for rng in _rngs():
                    with self.subTest(seeing=seeing, layers=name, rng=rng is not None):
                        h, w, info = generateProfile(
                            r0, theta0, gl if useGl else None, layers, rng)
                        self.assertAlmostEqual(
                            theta0FromProfile(h, w, r0) / theta0, 1, delta=RTOL)
                        if useGl:
                            self.assertAlmostEqual(_glFraction(h, w), gl, delta=RTOL)
                        self.assertAlmostEqual(info['glFraction'], _glFraction(h, w), delta=RTOL)

    def test_theta0_only_is_exact_and_info_matches_weights(self):
        """Mode A: theta0 exact; info['glFraction'] equals the GL fraction of the weights."""
        self._check(useGl=False)

    def test_theta0_and_gl_fraction_both_exact(self):
        """Mode B: theta0 and GL fraction both exact (includes the p90 jet-fallback case)."""
        self._check(useGl=True)

    def test_p90_condition_mode_b_median_succeeds(self):
        """Seeing 1.05, theta0 1.35, GL 0.291 (jet power hits 0, next bump takes over)."""
        r0 = seeingToR0(1.05)
        h, w, info = generateProfile(r0, 1.35, 0.291, DEFAULT_LAYERS)
        self.assertAlmostEqual(theta0FromProfile(h, w, r0) / 1.35, 1, delta=RTOL)
        self.assertAlmostEqual(_glFraction(h, w), 0.291, delta=RTOL)
        self.assertTrue(np.all(info['bumpPowers'] >= -1e-15))

    def test_gl_fraction_only_is_exact(self):
        """theta0=None: GL fraction exact."""
        for gl in (0.3, 0.6, 0.8):
            for name, layers in LAYER_SETS.items():
                for rng in _rngs():
                    with self.subTest(gl=gl, layers=name):
                        h, w, _ = generateProfile(0.15, None, gl, layers, rng)
                        self.assertAlmostEqual(_glFraction(h, w), gl, delta=RTOL)

    def test_no_constraints_gives_median_gl_fraction(self):
        """Neither theta0 nor GL fraction: GL fraction == F_GL_MEDIAN."""
        h, w, info = generateProfile(0.15)
        self.assertAlmostEqual(_glFraction(h, w), F_GL_MEDIAN, delta=RTOL)
        self.assertAlmostEqual(info['glFraction'], F_GL_MEDIAN, delta=RTOL)


class TestGenerateProfileWeights(unittest.TestCase):

    def test_weights_nonnegative_normalized_and_right_length(self):
        """Weights >= 0, sum to 1, one per requested layer."""
        for seeing, theta0, gl in CONDITIONS:
            r0 = seeingToR0(seeing)
            for name, layers in LAYER_SETS.items():
                n = len(GRID) if layers is None else len(layers)
                for useGl in (False, True):
                    for rng in _rngs():
                        with self.subTest(seeing=seeing, layers=name, useGl=useGl):
                            h, w, _ = generateProfile(
                                r0, theta0, gl if useGl else None, layers, rng)
                            self.assertEqual(len(w), n)
                            self.assertEqual(len(h), n)
                            self.assertGreaterEqual(w.min(), -1e-15)
                            self.assertAlmostEqual(w.sum(), 1, delta=1e-12)

    def test_returned_heights_are_requested_layers(self):
        """Returned heights equal the requested layers (GRID for None)."""
        r0 = seeingToR0(0.725)
        h, _, _ = generateProfile(r0, 2.0, layerHeights=None)
        np.testing.assert_array_equal(h, GRID)
        h, _, _ = generateProfile(r0, 2.0, layerHeights=LAYER_SETS['coarse'])
        np.testing.assert_array_equal(h, LAYER_SETS['coarse'])


class TestGenerateProfileErrors(unittest.TestCase):

    def test_gl_fraction_below_minimum_raises(self):
        """glFraction < F_GL_MIN is not reachable."""
        with self.assertRaises(ValueError):
            generateProfile(0.15, glFraction=0.01)
        with self.assertRaises(ValueError):
            generateProfile(0.15, glFraction=F_GL_MIN * 0.9)

    def test_gl_fraction_above_one_raises(self):
        """glFraction > 1 is not physical."""
        with self.assertRaises(ValueError):
            generateProfile(0.15, glFraction=1.01)

    def test_unreachable_theta0_raises(self):
        """theta0 far too large or far too small (median profile) raises, in modes A and B."""
        r0 = seeingToR0(0.725)
        for theta0 in (50., 0.05):
            for gl in (None, 0.6):
                with self.subTest(theta0=theta0, gl=gl):
                    with self.assertRaises(ValueError):
                        generateProfile(r0, theta0, gl, rng=None)

    def test_layer_matrix_needs_ground_and_high_layers(self):
        """layerMatrix raises if no layer <= 1000 m or none above."""
        with self.assertRaises(ValueError):
            layerMatrix([1500., 3000., 8000.])
        with self.assertRaises(ValueError):
            layerMatrix([0., 500., 1000.])


class TestGenerateProfileWarning(unittest.TestCase):

    def test_warns_when_gl_fraction_outside_observed_range(self):
        """GL fraction 0.97 (> F_GL_RANGE[1], reachable) emits a UserWarning."""
        self.assertGreater(0.97, F_GL_RANGE[1])
        with self.assertWarns(UserWarning):
            generateProfile(0.15, None, 0.97)

    def test_no_warning_inside_observed_range(self):
        """GL fraction 0.6 emits no warning."""
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            generateProfile(0.15, None, 0.6)


class TestGenerateProfileReproducibility(unittest.TestCase):

    def _draw(self, seed):
        return generateProfile(seeingToR0(0.725), 2.0, rng=np.random.default_rng(seed))[1]

    def test_same_seed_same_weights(self):
        """Same seed gives identical weights."""
        np.testing.assert_array_equal(self._draw(5), self._draw(5))

    def test_different_seeds_different_weights(self):
        """Different seeds give different weights."""
        self.assertFalse(np.array_equal(self._draw(5), self._draw(6)))

    def test_rng_none_is_deterministic(self):
        """rng=None gives the same median profile every time."""
        r0 = seeingToR0(0.725)
        w1 = generateProfile(r0, 2.0)[1]
        w2 = generateProfile(r0, 2.0)[1]
        np.testing.assert_array_equal(w1, w2)


class TestLogGrid(unittest.TestCase):

    def test_endpoints_length_and_monotonic(self):
        """First layer 0, last hTop, n layers, strictly increasing."""
        for n, h0, hTop in [(20, 1500., 20000.), (10, 1000., 25000.), (35, 500., 18000.)]:
            with self.subTest(n=n):
                g = logGrid(n, h0, hTop)
                self.assertEqual(len(g), n)
                self.assertEqual(g[0], 0)
                self.assertEqual(g[-1], hTop)
                self.assertTrue(np.all(np.diff(g) > 0))

    def test_spacing_increases_with_altitude(self):
        """Layer spacing grows with altitude (log spacing above h0)."""
        d = np.diff(logGrid())
        self.assertTrue(np.all(np.diff(d) >= 0))
        self.assertGreater(d[-1], 5 * d[0])


class TestLayerMatrix(unittest.TestCase):

    def test_one_per_column_and_power_conserved(self):
        """Each column has a single 1; M @ w sums to w.sum()."""
        for name, layers in LAYER_SETS.items():
            if layers is None:
                continue
            with self.subTest(layers=name):
                M = layerMatrix(layers)
                self.assertEqual(M.shape, (len(layers), len(GRID)))
                np.testing.assert_array_equal(M.sum(0), np.ones(len(GRID)))
                self.assertTrue(np.all((M == 0) | (M == 1)))
                w = np.random.default_rng(0).uniform(size=len(GRID))
                self.assertAlmostEqual((M @ w).sum(), w.sum(), delta=1e-12)

    def test_ground_layer_separation(self):
        """Bins <= H_GL only feed layers <= H_GL and bins above only layers above."""
        layers = np.asarray(LAYER_SETS['coarse'], dtype=float)
        M = layerMatrix(layers)
        low = GRID <= H_GL
        self.assertEqual(M[layers > H_GL][:, low].sum(), 0)
        self.assertEqual(M[layers <= H_GL][:, ~low].sum(), 0)


class TestLayerWind(unittest.TestCase):

    def setUp(self):
        self.layers = LAYER_SETS['coarse']
        _, self.w, _ = generateProfile(seeingToR0(0.725), 2.0, layerHeights=None)
        self.M = layerMatrix(self.layers)

    def test_uniform_speed_is_preserved(self):
        """Uniform wind speed v0 gives v0 on every layer."""
        v0 = 12.3
        speed, _ = layerWind(np.full(len(GRID), v0), np.zeros(len(GRID)), self.layers, self.w)
        np.testing.assert_allclose(speed, v0, rtol=RTOL)

    def test_uniform_direction_is_preserved(self):
        """Uniform wind direction d0 gives d0 on every layer."""
        d0 = 77.
        _, direction = layerWind(np.ones(len(GRID)), np.full(len(GRID), d0), self.layers, self.w)
        np.testing.assert_allclose(direction, d0, atol=1e-8)

    def test_direction_is_circular_mean(self):
        """350 and 10 deg with equal weights average to 0 deg, not 180."""
        layers = [0., 5000.]
        w = np.zeros(len(GRID))
        i, j = 0, 1                      # both bins below H_GL -> same layer
        w[[i, j]] = 0.5
        ang = np.full(len(GRID), 0.)
        ang[i], ang[j] = 350., 10.
        _, direction = layerWind(np.ones(len(GRID)), ang, layers, w + 1e-30)
        self.assertAlmostEqual((direction[0] + 180) % 360 - 180, 0, delta=1e-8)

    def test_tau0_preserved(self):
        """Cn2-weighted mean of v^(5/3) is the same on GRID and on the layers."""
        rng = np.random.default_rng(0)
        v = rng.uniform(2, 40, len(GRID))
        speed, _ = layerWind(v, np.zeros(len(GRID)), self.layers, self.w)
        wLayers = self.M @ self.w
        np.testing.assert_allclose(wLayers @ speed ** (5 / 3), self.w @ v ** (5 / 3), rtol=RTOL)


class TestSeeingR0(unittest.TestCase):

    def test_inverse(self):
        """seeingToR0 and r0ToSeeing are inverse of each other."""
        for s in (0.4, 0.8, 1.5):
            self.assertAlmostEqual(r0ToSeeing(seeingToR0(s)) / s, 1, delta=1e-12)

    def test_known_value(self):
        """seeing = 0.98 * 0.5e-6 / 0.1 rad corresponds to r0 = 0.1 m."""
        seeing = 0.98 * 0.5e-6 / 0.1 * 180 / np.pi * 3600
        self.assertAlmostEqual(seeingToR0(seeing), 0.1, delta=1e-12)


if __name__ == '__main__':
    unittest.main()
