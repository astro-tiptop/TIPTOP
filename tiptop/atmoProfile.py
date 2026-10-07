"""
Parametric Cn2 profiles with given integrated parameters, coherent with the
Stereo-SCIDAR Paranal statistics.

Model of a normalized profile w(h) (sum = 1), with f the GL fraction (Cn2 below 1 km)::

    w = GL bump      share g0 in the ground bin, the rest in a fixed tail up to 1 km
      + background   fixed shape, B (1 - f) / (1 - b_gl)   (its part above 1 km is B (1 - f))
      + 5 bumps      low [1-4 km], lowmid [3-7], mid [5-10], jet [10-17], high [17-25]:
                     Gaussians with a fixed FWHM per window, total power (1 - B)(1 - f)

g0, bump heights and bump power shares are drawn from the Stereo-SCIDAR statistics
(Gaussian copula), or set to their medians. theta0 fixes mu = sum(w h^5/3), because
theta0 = 0.314 r0 / mu^(3/5); the constraints are solved in closed form (mu is linear
in the free parameter)::

    theta0 only        -> f is solved;
    theta0 + f         -> the jet power is solved, exchanged with the other bumps; if the
                          jet reaches a limit, the next bump from the ground up takes over;
    f only / neither   -> theta0 is an output (f = Stereo-SCIDAR median if not given).

The model constants are in data/cn2_bump_model.npz, computed by bump_statistics.py
(TIPTOP_scripts/atmo) from the Stereo-SCIDAR Paranal database 2016-2021. All the
quantities are at zenith and 500 nm.
"""
import warnings
from pathlib import Path

import numpy as np
from scipy.stats import norm

RAD2ARCSEC = 180 / np.pi * 3600
WVL_REF = 500e-9          # reference wavelength of r0, seeing and theta0 [m]

# --- model constants
_m = np.load(Path(__file__).resolve().parent / 'data' / 'cn2_bump_model.npz')
GRID = _m['h']                                  # 250 m grid of the model [m]
H_GL = float(_m['h_gl'])                        # GL top [m]
BUMP_NAMES = [str(n) for n in _m['names']]
_FWHM = _m['fwhm']                              # bump FWHM per window [m]
_S_BG, _S_GL_TAIL = _m['s_bg'], _m['s_gl_tail']  # unit-sum shapes on GRID
_S_GROUND = (GRID == 0).astype(float)
_B, _B_GL = float(_m['B']), float(_m['b_gl'])
F_GL_MEDIAN = float(_m['f_gl_median'])
F_GL_RANGE = tuple(_m['f_gl_range'])            # p1-p99 observed at Paranal
_Q_PROBS, _Q_G0, _Q_H, _Q_P = _m['q_probs'], _m['q_g0'], _m['q_h'], _m['q_p']
_R_H, _R_P = _m['R_h'], _m['R_p']
_N_BUMPS = len(BUMP_NAMES)
_JET = BUMP_NAMES.index('jet')
# smallest f: below it the GL bump would be negative (background alone exceeds f)
F_GL_MIN = _B * _B_GL / (1 - _B_GL + _B * _B_GL)
# theta0 + f: the jet moves first, then the other bumps from the ground up
_KNOB_ORDER = [_JET] + [k for k in range(_N_BUMPS) if k != _JET]


def logGrid(n: int = 20, h0: float = 1500., hTop: float = 20000.) -> np.ndarray:
    """Layer altitudes uniform in log(1 + h/h0), from 0 to hTop.

    Spacing ~ (h + h0): about linear below h0, logarithmic above. With the defaults,
    ~250 m near the ground and ~2.5 km near 20 km. On MORFEO, 20 such layers give the
    same Strehl ratio as 40 uniform layers within 0.004 on average (TIPTOP_scripts/atmo/
    compare_layers.py, 5 conditions of Agapito+ 2021).

    Parameters
    ----------
    n : int
        Number of layers.
    h0 : float
        Transition altitude between linear and logarithmic spacing [m].
    hTop : float
        Altitude of the highest layer [m].

    Returns
    -------
    numpy.ndarray
        Layer altitudes [m], rounded to 10 m.
    """
    u = np.linspace(0, np.log(1 + hTop / h0), n)
    return np.round(h0 * np.expm1(u), -1)


DEFAULT_LAYERS = logGrid()


def layerMatrix(layerHeights) -> np.ndarray:
    """Matrix moving each 250 m bin of GRID to the nearest layer.

    Bins below (above) H_GL only go to layers below (above) H_GL, so that the GL
    fraction is the same on GRID and on the layers.

    Parameters
    ----------
    layerHeights : array_like
        Layer altitudes [m], at least one <= H_GL and one above.

    Returns
    -------
    numpy.ndarray
        Shape (n_layers, len(GRID)), entries 0 or 1, one 1 per column.
    """
    lh = np.asarray(layerHeights, dtype=float)
    side = {True: np.flatnonzero(lh <= H_GL), False: np.flatnonzero(lh > H_GL)}
    if len(side[True]) == 0 or len(side[False]) == 0:
        raise ValueError(f'layerHeights needs at least one layer <= {H_GL:.0f} m and one above')
    M = np.zeros((len(lh), len(GRID)))
    for i, hb in enumerate(GRID):
        cand = side[bool(hb <= H_GL)]
        M[cand[np.argmin(np.abs(lh[cand] - hb))], i] = 1.
    return M


def seeingToR0(seeing: float, wvl: float = WVL_REF) -> float:
    """Fried parameter [m] from the seeing FWHM [arcsec]: seeing = 0.98 wvl / r0."""
    return 0.98 * wvl / (seeing / RAD2ARCSEC)


def r0ToSeeing(r0: float, wvl: float = WVL_REF) -> float:
    """Seeing FWHM [arcsec] from the Fried parameter [m]: seeing = 0.98 wvl / r0."""
    return 0.98 * wvl / r0 * RAD2ARCSEC


def theta0FromProfile(heights, weights, r0: float) -> float:
    """Isoplanatic angle [arcsec]: theta0 = 0.314 r0 / (sum w h^5/3)^(3/5).

    Parameters
    ----------
    heights : array_like
        Layer altitudes [m].
    weights : array_like
        Cn2 weights, normalized to sum 1.
    r0 : float
        Fried parameter at the same wavelength [m].
    """
    h, w = np.asarray(heights, dtype=float), np.asarray(weights, dtype=float)
    return 0.314 * r0 / (w @ h ** (5 / 3)) ** (3 / 5) * RAD2ARCSEC


def generateProfile(r0: float, theta0: float = None, glFraction: float = None,
                    layerHeights=DEFAULT_LAYERS, rng: np.random.Generator = None,
                    maxTries: int = 100):
    """Generate a normalized Cn2 profile with given theta0 and/or GL fraction.

    See the module docstring for the model.

    Parameters
    ----------
    r0 : float
        Fried parameter at 500 nm, zenith [m]. Only used to convert theta0.
    theta0 : float, optional
        Isoplanatic angle at 500 nm, zenith [arcsec].
    glFraction : float, optional
        Fraction of the Cn2 below 1 km.
    layerHeights : array_like or None
        Output layer altitudes [m] (at least one <= 1 km and one above); default: the
        20 layers of logGrid(). theta0 and glFraction are imposed on these layers.
        None returns the 250 m grid GRID.
    rng : numpy.random.Generator, optional
        Random draw of g0, bump heights and power shares. None gives the median profile.
    maxTries : int
        Random draws incompatible with (theta0, glFraction) are redrawn up to this
        many times.

    Returns
    -------
    heights : numpy.ndarray
        Layer altitudes [m].
    weights : numpy.ndarray
        Normalized Cn2 weights (sum = 1).
    info : dict
        'glFraction', 'g0', 'bumpHeights' [m] and 'bumpPowers'.

    Raises
    ------
    ValueError
        If (theta0, glFraction) cannot be reached by the model.

    Warns
    -----
    UserWarning
        If the GL fraction is outside the range observed at Paranal (p1-p99).
    """
    if glFraction is not None and not F_GL_MIN <= glFraction <= 1:
        raise ValueError(f'glFraction = {glFraction:.3f} outside [{F_GL_MIN:.3f}, 1]: below '
                         f'{F_GL_MIN:.3f} the background alone puts more Cn2 below 1 km')
    muTarget = None if theta0 is None else (0.314 * r0 / (theta0 / RAD2ARCSEC)) ** (5 / 3)
    h = GRID if layerHeights is None else np.asarray(layerHeights, dtype=float)
    toLayers = np.eye(len(GRID)) if layerHeights is None else layerMatrix(h)
    h53 = h ** (5 / 3)
    sBg = toLayers @ _S_BG

    for _ in range(maxTries):
        g0, bumpHeights, shares = _drawBumps(rng)
        sGl = toLayers @ (g0 * _S_GROUND + (1 - g0) * _S_GL_TAIL)
        bumps = _bumpShapes(bumpHeights) @ toLayers.T
        f, powers, problem = _solvePowers(muTarget, glFraction, shares,
                                          sGl @ h53, sBg @ h53, bumps @ h53)
        if problem is None or rng is None:   # the median profile cannot be redrawn
            break
    if problem is not None:
        raise ValueError(f'r0 = {r0:.3f} m, theta0 = {theta0}", glFraction = {glFraction}: '
                         f'{problem}')
    if not F_GL_RANGE[0] <= f <= F_GL_RANGE[1]:
        warnings.warn(f'GL fraction {f:.3f} outside the range observed at Paranal '
                      f'[{F_GL_RANGE[0]:.2f}, {F_GL_RANGE[1]:.2f}] (p1-p99)')

    pGl, pBg, _ = _components(f, shares)
    weights = pGl * sGl + pBg * sBg + powers @ bumps
    return h, weights, {'glFraction': f, 'g0': g0, 'bumpHeights': bumpHeights,
                        'bumpPowers': powers}


def layerWind(windSpeed, windDirection, layerHeights, weights) -> tuple:
    """Wind speed and direction on the layers, averaged with Cn2 weights.

    Speed: (sum w v^5/3 / sum w)^(3/5) in each layer, so that tau0 is preserved.
    Direction: Cn2-weighted circular mean.

    Parameters
    ----------
    windSpeed, windDirection : array_like
        Wind speed [m/s] and direction [deg] on GRID (250 m).
    layerHeights : array_like
        Layer altitudes [m].
    weights : array_like
        Cn2 weights on GRID (e.g. generateProfile(..., layerHeights=None)).

    Returns
    -------
    speed, direction : numpy.ndarray
        Wind speed [m/s] and direction [deg] on the layers.
    """
    w = np.asarray(weights, dtype=float)
    v = np.asarray(windSpeed, dtype=float)
    ang = np.deg2rad(np.asarray(windDirection, dtype=float))
    M = layerMatrix(layerHeights)
    speed = ((M @ (w * v ** (5 / 3))) / (M @ w)) ** (3 / 5)
    direction = np.rad2deg(np.arctan2(M @ (w * np.sin(ang)), M @ (w * np.cos(ang))))
    return speed, direction


def _drawBumps(rng):
    """g0, bump heights [m] and power shares (sum 1); medians if rng is None."""
    if rng is None:
        uG, uH, uP = 0.5, np.full(_N_BUMPS, 0.5), np.full(_N_BUMPS, 0.5)
    else:
        uG = rng.uniform()
        uH = norm.cdf(rng.multivariate_normal(np.zeros(_N_BUMPS), _R_H))
        uP = norm.cdf(rng.multivariate_normal(np.zeros(_N_BUMPS), _R_P))
    g0 = np.interp(uG, _Q_PROBS, _Q_G0)
    heights = np.array([np.interp(uH[k], _Q_PROBS, _Q_H[:, k]) for k in range(_N_BUMPS)])
    powers = np.array([np.interp(uP[k], _Q_PROBS, _Q_P[:, k]) for k in range(_N_BUMPS)])
    return g0, heights, powers / powers.sum()


def _bumpShapes(heights):
    """Gaussians with the FWHM of their window, zero below H_GL, unit sum on GRID."""
    sigma = _FWHM / (2 * np.sqrt(2 * np.log(2)))
    g = np.exp(-0.5 * ((GRID - heights[:, None]) / sigma[:, None]) ** 2) * (GRID > H_GL)
    return g / g.sum(1, keepdims=True)


def _components(f, shares):
    """Powers of GL bump, background and bumps for a given GL fraction f."""
    pBg = _B * (1 - f) / (1 - _B_GL)
    return f - _B_GL * pBg, pBg, (1 - _B) * (1 - f) * shares


def _solvePowers(muTarget, glFraction, shares, mGl, mBg, mBumps):
    """Apply the theta0 / GL fraction constraints.

    mGl, mBg, mBumps: contribution to mu of the unit GL bump, background and bumps.
    Returns f, bump powers and a problem description (None if solved).
    """
    def mu(f):
        pGl, pBg, pK = _components(f, shares)
        return pGl * mGl + pBg * mBg + pK @ mBumps

    if muTarget is None:
        f = F_GL_MEDIAN if glFraction is None else glFraction
        return f, _components(f, shares)[2], None

    if glFraction is None:
        # mu(f) = mu(0) + f (mu(1) - mu(0))
        f = (muTarget - mu(0)) / (mu(1) - mu(0))
        if not F_GL_MIN <= f <= 1:
            return f, None, f'needs a GL fraction {f:.3f}, outside [{F_GL_MIN:.3f}, 1]'
        return f, _components(f, shares)[2], None

    # f fixed: knob k moves power x between bump k and the bumps not fixed yet (which
    # keep their relative shares); mu is linear in x. If x falls outside [0, available
    # power], bump k is fixed at that limit and the next knob goes on.
    pGl, pBg, pK = _components(glFraction, shares)
    powers = np.full(_N_BUMPS, np.nan)              # nan = not fixed yet
    pFree = pK.sum()
    muFree = muTarget - pGl * mGl - pBg * mBg
    for k in _KNOB_ORDER:
        others = np.isnan(powers)
        others[k] = False
        if not others.any():
            break
        q = np.where(others, shares, 0.) / shares[others].sum()
        mOthers = q @ mBumps
        x = (muFree - pFree * mOthers) / (mBumps[k] - mOthers)
        if 0 <= x <= pFree:
            return glFraction, np.where(others, (pFree - x) * q, np.nan_to_num(powers)) \
                + np.eye(_N_BUMPS)[k] * x, None
        powers[k] = np.clip(x, 0, pFree)
        pFree -= powers[k]
        muFree -= powers[k] * mBumps[k]
    return glFraction, None, f'theta0 not reachable with GL fraction {glFraction:.3f}'
