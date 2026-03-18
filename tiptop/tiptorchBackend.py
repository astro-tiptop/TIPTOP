"""
tiptorchBackend.py
==================
Adapter layer that lets TIPTOP's baseSimulation use TipTorch as an alternative
HO-PSD computation engine instead of P3's fourierModel.

Public API
----------
build_tiptorch_config(my_data_map,
                      zenith_science, azimuth_science,
                      zenith_ngs, azimuth_ngs,
                      device, dtype)
    -> (config, ao_type, n_sources)

build_freq_proxy_from_tiptorch(model)
    -> TipTorchFreqProxy
        A lightweight namespace that mimics the P3 frequencyDomain attributes
        consumed by baseSimulation after the HO PSD has been produced.

run_tiptorch_backend(my_data_map, LOisOn,
                     zenith_science, azimuth_science,
                     zenith_ngs, azimuth_ngs,
                     device, dtype, tiptorch_kwargs)
    -> (psd_numpy, model, freq_proxy)
        psd_numpy  : numpy array (N, N, N_src) in nm^2, P3 convention
        model      : live TipTorch instance (for inspection / fitting)
        freq_proxy : TipTorchFreqProxy (replaces P3's frequencyDomain)

Design notes
------------
When TipTorch is the backend, baseSimulation no longer calls
fourierModel.initComputations().  Instead it:
  1. Instantiates fourierModel(doComputations=False) to obtain fao.ao
     (telescope pupil, P3 atmosphere object for open-loop PSF, etc.).
  2. Calls run_tiptorch_backend() to compute the HO PSD.
  3. Assigns fao.freq = freq_proxy  and  fao.PSD = psd_numpy.
  4. Sets fao.ao.atm.wvl = freq_proxy.wvlRef so computeOL_PSD() works.

All downstream TIPTOP stages (MASTSEL, LO convolution, OL/DL PSFs, file
output) remain completely unchanged.
"""

import numpy as np
import torch
import warnings


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _to_tensor(value, device, dtype):
    arr = np.atleast_1d(np.asarray(value, dtype=np.float64))
    return torch.as_tensor(arr, device=device, dtype=dtype)


def _first(value):
    """Return value[0] if list/tuple, else value itself."""
    return value[0] if isinstance(value, (list, tuple)) else value


def _as_list(value):
    return list(value) if isinstance(value, (list, tuple)) else [value]


# ---------------------------------------------------------------------------
# AO-type detection
# ---------------------------------------------------------------------------

def detect_ao_type(my_data_map):
    """
    Infer the TipTorch AO_type string from the TIPTOP parameter map.
    Returns one of: 'SCAO', 'LTAO', 'MCAO', 'GLAO'.
    """
    has_lgs = 'sources_LGS' in my_data_map
    n_dm = len(_as_list(my_data_map.get('DM', {}).get('DmHeights', [0.0])))
    if has_lgs:
        return 'MCAO' if n_dm > 1 else 'LTAO'
    return 'SCAO'


# ---------------------------------------------------------------------------
# TipTorch config builder
# ---------------------------------------------------------------------------

def build_tiptorch_config(my_data_map,
                          zenith_science,
                          azimuth_science,
                          zenith_ngs=None,
                          azimuth_ngs=None,
                          device=torch.device('cpu'),
                          dtype=torch.float32):
    """
    Convert a TIPTOP my_data_map dict into the tensor-based config dict that
    TipTorch.__init__() expects.

    Science and (optional) NGS directions are merged into a single source list
    so that one TipTorch forward pass produces PSDs for all directions.

    Parameters
    ----------
    my_data_map     : dict
    zenith_science  : array-like, arcsec
    azimuth_science : array-like, deg
    zenith_ngs      : array-like or None, arcsec
    azimuth_ngs     : array-like or None, deg
    device, dtype   : torch

    Returns
    -------
    config    : dict  (ready for TipTorch.__init__)
    ao_type   : str
    n_sources : int   (science + ngs)
    """
    tel    = my_data_map.get('telescope', {})
    atm    = my_data_map.get('atmosphere', {})
    dm     = my_data_map.get('DM', {})
    s_sci  = my_data_map.get('sources_science', {})
    ss     = my_data_map.get('sensor_science', {})
    rtc    = my_data_map.get('RTC', {})

    # HO WFS / guide-star section (multiple possible key names across configs)
    s_ho   = (my_data_map.get('sources_HO')
              or my_data_map.get('sources_LGS')
              or my_data_map.get('sources_NGS') or {})
    sen_ho = (my_data_map.get('sensor_HO')
              or my_data_map.get('sensor_WFS') or {})

    ao_type = detect_ao_type(my_data_map)

    # ---- Combined source list -----------------------------------------------
    zen_sci = _as_list(zenith_science)
    az_sci  = _as_list(azimuth_science)
    if zenith_ngs is not None:
        zen_all = zen_sci + _as_list(zenith_ngs)
        az_all  = az_sci  + _as_list(azimuth_ngs)
    else:
        zen_all, az_all = zen_sci, az_sci
    n_sources = len(zen_all)

    t = lambda v: _to_tensor(v, device, dtype)

    # ---- Science wavelength(s) ---------------------------------------------
    wvl_sci = _as_list(s_sci.get('Wavelength', [0.5e-6]))

    # ---- Atmosphere --------------------------------------------------------
    atm_wvl = float(_first(atm.get('Wavelength', 500e-9)))

    # ---- HO guide stars ----------------------------------------------------
    ho_zen = _as_list(s_ho.get('Zenith',  [0.0]))
    ho_az  = _as_list(s_ho.get('Azimuth', [0.0]))
    ho_h   = _as_list(s_ho.get('Height',  [0.0]))
    ho_wvl = float(_first(s_ho.get('Wavelength', 589e-9)))

    # ---- DM ----------------------------------------------------------------
    dm_pitches = _as_list(dm.get('DmPitchs',           [0.5]))
    dm_heights = _as_list(dm.get('DmHeights',           [0.0]))
    dm_opt_zen = _as_list(dm.get('OptimizationZenith',  [0.0]))
    dm_opt_az  = _as_list(dm.get('OptimizationAzimuth', [0.0]))
    dm_opt_wt  = _as_list(dm.get('OptimizationWeight',  [1.0]))

    # ---- HO sensor ---------------------------------------------------------
    D       = float(tel['TelescopeDiameter'])
    d_sub   = float(_first(sen_ho.get('SizeLenslets',      dm_pitches[0])))
    n_sub   = int(_first(sen_ho.get('NumberLenslets',
                         int(D / dm_pitches[0]))))
    ho_rate = float(_first(rtc.get('SensorFrameRate_HO',   500.0)))
    wfs_clk = float(_first(sen_ho.get('ClockRate',          ho_rate)))
    wfs_fov = float(_first(sen_ho.get('FieldOfView',        1.0)))
    wfs_ron = float(_first(sen_ho.get('SigmaRON',           0.5)))
    wfs_exc = float(_first(sen_ho.get('ExcessNoiseFactor',  1.0)))
    wfs_nph = float(_first(sen_ho.get('NumberPhotons',    1000.0)))
    wfs_ps  = float(_first(sen_ho.get('PixelScale',        500.0)))
    wfs_fwhm = _as_list(sen_ho.get('SpotFWHM', [[0.0, 0.0, 0.0]]))
    if not isinstance(wfs_fwhm[0], (list, tuple)):
        wfs_fwhm = [[wfs_fwhm[0], wfs_fwhm[0], 0.0]]
    wfs_alg = str(sen_ho.get('Algorithm',         'wcog')).lower()
    wfs_wcg = float(_first(sen_ho.get('WindowRadiusWCoG',  2.0)))
    wfs_typ = str(sen_ho.get('WfsType',   'Shack-Hartmann'))

    # ---- RTC ---------------------------------------------------------------
    ho_delay = float(_first(rtc.get('LoopDelaySteps_HO',  2.0)))
    ho_gain  = float(_first(rtc.get('LoopGain_HO',        0.3)))

    # ---- Telescope ---------------------------------------------------------
    zenith_ang = float(_first(tel.get('ZenithAngle',   0.0)))
    pupil_path = tel.get('PathPupil',    None)
    apod_path  = tel.get('PathApodizer', None)
    pupil_ang  = float(_first(tel.get('PupilAngle',    0.0)))

    # ---- Science sensor ----------------------------------------------------
    ps_mas = float(_first(ss.get('PixelScale',  1.0)))
    fov_px = int(_first(ss.get('FieldOfView', 128)))

    config = {
        'NumberSources': n_sources,

        'telescope': {
            'TelescopeDiameter': D,
            'ZenithAngle':  t([zenith_ang]),
            'PathPupil':    pupil_path,
            'PathApodizer': apod_path,
            'PupilAngle':   pupil_ang,
        },

        'atmosphere': {
            'Wavelength':    atm_wvl,
            'Seeing':        t([float(_first(atm.get('Seeing', 0.65)))]),
            'L0':            t(_as_list(atm.get('OuterScale',
                                                atm.get('L0', [25.0])))),
            'Cn2Weights':    t(_as_list(atm.get('Cn2Weights', [1.0]))),
            'Cn2Heights':    t(_as_list(atm.get('Cn2Heights', [0.0]))),
            'WindSpeed':     t(_as_list(atm.get('WindSpeed',  [10.0]))),
            'WindDirection': t(_as_list(atm.get('WindDirection', [0.0]))),
        },

        'sources_science': {
            'Wavelength': t(wvl_sci),   # [N_wvl]
            'Zenith':     t(zen_all),   # [N_src] arcsec
            'Azimuth':    t(az_all),    # [N_src] deg
        },

        'sources_HO': {
            'Wavelength': t([ho_wvl]),
            'Height':     t(ho_h),
            'Zenith':     t(ho_zen),
            'Azimuth':    t(ho_az),
        },

        'sensor_science': {
            'PixelScale':  ps_mas,
            'FieldOfView': fov_px,
        },

        'DM': {
            'DmPitchs':            t(dm_pitches),
            'DmHeights':           t(dm_heights),
            'OptimizationZenith':  t(dm_opt_zen),
            'OptimizationAzimuth': t(dm_opt_az),
            'OptimizationWeight':  t(dm_opt_wt),
        },

        'sensor_HO': {
            'SizeLenslets':      t([d_sub]),
            'NumberLenslets':    t([n_sub]),
            'ClockRate':         t([wfs_clk]),
            'FieldOfView':       wfs_fov,
            'SigmaRON':          t([wfs_ron]),
            'PixelScale':        t([wfs_ps]),
            'SpotFWHM':          t(wfs_fwhm),
            'ExcessNoiseFactor': t([wfs_exc]),
            'NumberPhotons':     t([wfs_nph]),
            'Algorithm':         wfs_alg,
            'WindowRadiusWCoG':  t([wfs_wcg]),
            'WfsType':           wfs_typ,
        },

        'RTC': {
            'SensorFrameRate_HO': t([ho_rate]),
            'LoopDelaySteps_HO':  t([ho_delay]),
            'LoopGain_HO':        t([ho_gain]),
        },
    }

    return config, ao_type, n_sources


# ---------------------------------------------------------------------------
# Frequency-domain proxy
# ---------------------------------------------------------------------------

class TipTorchFreqProxy:
    """
    Lightweight namespace that stands in for P3's ``frequencyDomain`` object
    inside baseSimulation when TipTorch is the HO-PSD backend.

    Only the attributes actually consumed by baseSimulation are populated.

    Attribute mapping (P3 frequencyDomain  <->  TipTorch InitGrids)
    ──────────────────────────────────────────────────────────────────
    P3 frequencyDomain              TipTorch model
    ────────────────────────────    ────────────────────────────────────
    kRef_   = ceil(2/samp)          nOtf / N_pix  (eff. oversampling)
    PSDstep = 1/(D * samp * k_)     2*kc / nOtf_AO  (= dk post-correction)
    kcMax_  = max(1/(2*pitch))      kc.max().item()
    resAO   = int(2*kcMax_/step)    nOtf_AO
    wvlRef  = min(science wvl)      wvl_atm  (PSD normalised at atm wvl)
    psInMas = ao.cam.psInMas        config['sensor_science']['PixelScale']
    k2_     = freq_array(nOtf)**2   reconstructed from nOtf and PSDstep

    Derivation that PSDstep = 2*kc/nOtf_AO
    ─────────────────────────────────────────
    P3:       PSDstep = psInMas / (wvl * rad2mas * k_)
                      = 1 / (D * sampling)   where sampling = k_ * samp
    TipTorch: dk_init = 1 / D / sampling_min
    Mask gives nOtf_AO = round(2*kc / dk_init)
    => dk_corrected = 2*kc / nOtf_AO = dk_init = PSDstep  ✓
    """

    def __repr__(self):
        attrs = ['kRef_', 'PSDstep', 'kcMax_', 'resAO', 'wvlRef', 'psInMas']
        lines = ['TipTorchFreqProxy(']
        for a in attrs:
            if hasattr(self, a):
                lines.append(f'  {a} = {getattr(self, a)}')
        lines.append(')')
        return '\n'.join(lines)


def build_freq_proxy_from_tiptorch(model):
    """
    Populate a TipTorchFreqProxy from an initialised TipTorch model.

    Must be called *after* TipTorch.__init__() has completed (i.e. after
    InitGrids has run), so that nOtf, nOtf_AO, dk, kc, etc. are available.

    Parameters
    ----------
    model : TipTorch
        Fully initialised TipTorch instance.

    Returns
    -------
    proxy : TipTorchFreqProxy
    """
    proxy = TipTorchFreqProxy()

    N    = int(model.nOtf)
    Npix = int(model.N_pix)

    # ------------------------------------------------------------------
    # 1.  Cut-off frequency and AO-corrected region
    #     kc may be a multi-element tensor (multiple DMs); take the max.
    # ------------------------------------------------------------------
    kc_val  = float(model.kc.max().item()) if model.kc.numel() > 1 \
              else float(model.kc.item())
    nOtf_AO = int(model.nOtf_AO)

    proxy.kcMax_ = kc_val
    proxy.resAO  = nOtf_AO

    # ------------------------------------------------------------------
    # 2.  Spatial frequency step  (see class docstring for derivation)
    # ------------------------------------------------------------------
    proxy.PSDstep = kc_val * 2.0 / nOtf_AO

    # ------------------------------------------------------------------
    # 3.  Oversampling factor
    #     kRef_ = nOtf / N_pix.  Due to TipTorch's _to_odd() this may
    #     not be an exact integer; baseSimulation casts it with int().
    # ------------------------------------------------------------------
    proxy.kRef_ = N / Npix

    # ------------------------------------------------------------------
    # 4.  Reference wavelength
    #     TipTorch normalises its PSD at wvl_atm.  The caller must also
    #     set fao.ao.atm.wvl = proxy.wvlRef so computeOL_PSD() uses the
    #     same wavelength when evaluating the P3 atmospheric spectrum.
    # ------------------------------------------------------------------
    proxy.wvlRef = float(model.wvl_atm)

    # ------------------------------------------------------------------
    # 5.  Pixel scale  — stored as a 1-element list so psInMas[0] works
    # ------------------------------------------------------------------
    proxy.psInMas = [float(model.psInMas)]

    # ------------------------------------------------------------------
    # 6.  Full k^2 spatial-frequency grid  (N x N, float64, m^{-2})
    #
    #     Replicates P3's freq_array(nOtf, offset=1e-10, L=PSDstep):
    #       np.mgrid[-N//2 : N//2 : N*1j] * PSDstep + 1e-10
    #     For odd N (guaranteed by TipTorch's _to_odd) this equals:
    #       (arange(N) - N//2) * PSDstep + 1e-10
    #
    #     Used by ngsPSF() and computeOL_PSD() as:
    #       k = np.sqrt(self.fao.freq.k2_)
    #       pf = FourierUtils.pistonFilter(diameter, k)
    # ------------------------------------------------------------------
    step   = proxy.PSDstep
    coords = (np.arange(N, dtype=np.float64) - N // 2) * step + 1e-10
    proxy.k2_ = (coords[:, None]**2 + coords[None, :]**2).astype(np.float64)

    return proxy


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_tiptorch_backend(my_data_map,
                         LOisOn,
                         zenith_science,
                         azimuth_science,
                         nNaturalGS_field=0,
                         zenith_ngs=None,
                         azimuth_ngs=None,
                         device=torch.device('cpu'),
                         dtype=torch.float32,
                         tiptorch_kwargs=None):
    """
    Run TipTorch to compute the HO PSD for all science + NGS directions.

    No P3 objects are needed by this function.

    Parameters
    ----------
    my_data_map      : dict   TIPTOP parameter dict
    LOisOn           : bool
    zenith_science   : array-like, arcsec
    azimuth_science  : array-like, deg
    nNaturalGS_field : int    number of NGS guide stars (0 when LOisOn=False)
    zenith_ngs       : array-like or None, arcsec
    azimuth_ngs      : array-like or None, deg
    device           : torch.device
    dtype            : torch.dtype
    tiptorch_kwargs  : dict or None   extra kwargs for TipTorch.__init__()

    Returns
    -------
    psd_p3     : numpy.ndarray, shape (N, N, N_src), nm^2
                 P3 array convention — transpose in baseSimulation gives
                 (N_src, N, N).
    model      : TipTorch   live model (for inspection / fitting / GPU)
    freq_proxy : TipTorchFreqProxy   drop-in replacement for fao.freq
    """
    try:
        from TipTorch import TipTorch
    except ImportError:
        raise ImportError(
            "The 'TipTorch' module could not be imported.  "
            "Ensure it is on sys.path before using backend='tiptorch'."
        )

    if tiptorch_kwargs is None:
        tiptorch_kwargs = {}

    # ------------------------------------------------------------------
    # 1.  Build config and instantiate TipTorch
    #     __init__ calls Update -> InitValues + InitGrids + InitPupils
    # ------------------------------------------------------------------
    zen_ngs = zenith_ngs  if (LOisOn and zenith_ngs  is not None) else None
    az_ngs  = azimuth_ngs if (LOisOn and azimuth_ngs is not None) else None

    config, ao_type, _ = build_tiptorch_config(
        my_data_map,
        zenith_science=zenith_science,
        azimuth_science=azimuth_science,
        zenith_ngs=zen_ngs,
        azimuth_ngs=az_ngs,
        device=device,
        dtype=dtype,
    )

    model = TipTorch(
        AO_config=config,
        AO_type=ao_type,
        device=device,
        dtype=dtype,
        **tiptorch_kwargs
    )

    # ------------------------------------------------------------------
    # 2.  Build the frequencyDomain proxy from the now-initialised grid
    #     (nOtf, nOtf_AO, kc, dk, etc. are set by InitGrids inside __init__)
    # ------------------------------------------------------------------
    freq_proxy = build_freq_proxy_from_tiptorch(model)

    # ------------------------------------------------------------------
    # 3.  Compute PSD  (N_src, N_wvl, nOtf, nOtf)  in nm^2
    # ------------------------------------------------------------------
    with torch.no_grad():
        psd_torch = model.ComputePSD()

    # ------------------------------------------------------------------
    # 4.  Convert to P3 array convention: (N, N, N_src)
    #     Only the first wavelength slice is taken — the PSD is at
    #     wvl_atm which matches P3's wvlRef convention.
    # ------------------------------------------------------------------
    psd_np = psd_torch.detach().cpu().numpy()   # (N_src, N_wvl, N, N)
    pupil_np = model.pupil.detach().cpu().numpy()   # (N, N)
    
    psd_np = psd_np[:, 0, :, :]                 # (N_src, N, N)
    psd_p3 = psd_np.transpose(1, 2, 0).astype(np.float64)  # (N, N, N_src)
    
    return psd_p3[:-1, :-1, :], model, freq_proxy, pupil_np