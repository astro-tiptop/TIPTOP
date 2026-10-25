import types
import numpy as np
import pytest

from mastsel import convolve, residualToSpectrum, Field
from tiptop.abstractSimulation import AbstractSimulation

N = 128
PS_MAS = 5.0
WVL = 1.65e-6
SIGMA0_PIX = 4.0


def _gaussian_field(sigma_pix: float) -> Field:
    """Unit-sum circular Gaussian PSF centered on the MASTSEL grid center."""
    probe = Field(WVL, N, 1.0, 'm')
    xp = probe.xp
    yy, xx = np.mgrid[:N, :N]
    c = N // 2
    g = np.exp(-((xx - c)**2 + (yy - c)**2) / (2 * sigma_pix**2))
    probe.sampling = xp.asarray(g / g.sum(), dtype=probe.sampling.dtype)
    return probe


def _kernel(theta: float, sx_pix: float, sy_pix: float) -> Field:
    """Kernel from an ellipse in mas, as built in baseSimulation._apply_final_convolution."""
    ellp = np.array([theta, sx_pix * PS_MAS, sy_pix * PS_MAS], dtype=np.float64)
    return residualToSpectrum(ellp, WVL, N, 1 / (N * PS_MAS))


def _host(f: Field) -> np.ndarray:
    s = f.sampling
    return s.get() if hasattr(s, 'get') else np.asarray(s)


def test_single_convolution_with_kernel_product_equals_two_convolutions():
    psf = _gaussian_field(SIGMA0_PIX)
    k1 = _kernel(0.4, 2.0, 1.0)
    k2 = _kernel(0.0, 1.5, 1.5)

    two = _host(convolve(convolve(psf, k1), k2))

    k12 = _kernel(0.4, 2.0, 1.0)
    k12.sampling = k12.sampling * k2.sampling
    one = _host(convolve(psf, k12))

    np.testing.assert_allclose(one, two, rtol=0, atol=1e-12 * two.max())


def test_second_moments_add_and_flux_is_one():
    psf = _gaussian_field(SIGMA0_PIX)
    s1x, s1y, s2x, s2y = 2.0, 3.0, 1.5, 2.5
    k = _kernel(0.0, s1x, s1y)
    k2 = _kernel(0.0, s2x, s2y)
    k.sampling = k.sampling * k2.sampling
    out = _host(convolve(psf, k))

    assert out.sum() == pytest.approx(1.0, abs=1e-9)
    yy, xx = np.mgrid[:N, :N]
    mx, my = (out * xx).sum(), (out * yy).sum()
    varx = (out * (xx - mx)**2).sum()
    vary = (out * (yy - my)**2).sum()
    assert varx == pytest.approx(SIGMA0_PIX**2 + s1x**2 + s2x**2, rel=1e-3)
    assert vary == pytest.approx(SIGMA0_PIX**2 + s1y**2 + s2y**2, rel=1e-3)


def _stub(cube: np.ndarray) -> types.SimpleNamespace:
    return types.SimpleNamespace(cubeResultsArray=cube, nWvl=1, psInMas=PS_MAS,
                                 SupSamp=None, verbose=False)


def _unit_gaussian_psf(sigma_pix: float = 6.0) -> np.ndarray:
    yy, xx = np.mgrid[:N, :N]
    g = np.exp(-((xx - N // 2)**2 + (yy - N // 2)**2) / (2 * sigma_pix**2))
    return g / g.sum()


def test_psf1d_independent_of_scale_and_input_unmodified():
    psf = _unit_gaussian_psf()
    cube_unit = psf[None].copy()
    cube_scaled = 7.0 * psf[None]
    backup = cube_scaled.copy()

    s_unit, s_scaled = _stub(cube_unit), _stub(cube_scaled)
    AbstractSimulation.computePSF1D(s_unit)
    AbstractSimulation.computePSF1D(s_scaled)

    np.testing.assert_allclose(s_scaled.psf1d, s_unit.psf1d, rtol=0, atol=1e-14)
    np.testing.assert_array_equal(cube_scaled, backup)
    np.testing.assert_array_equal(cube_unit, psf[None])


def test_psf1d_gaussian_profile_decreasing_and_peak_value():
    psf = _unit_gaussian_psf()
    s = _stub(psf[None].copy())
    AbstractSimulation.computePSF1D(s)
    prof = np.asarray(s.psf1d)[0, 0]

    assert np.all(np.diff(prof[:40]) < 0)
    # first bin averages the center pixel with its neighbours: within a few percent of the peak
    assert prof[0] == pytest.approx(psf.max(), rel=5e-2)
