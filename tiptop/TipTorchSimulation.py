"""
TIPTOP simulation backed by TipTorch.

TipTorch replaces P3 for all PSDs and for the science PSFs. MASTSEL still provides the low-order residual covariances
(MavisLO) and the NGS / Focus sensor spot PSFs. The public interface mirrors `tiptop.baseSimulation.baseSimulation`
so that `tiptop.py` and `asterismSimulation.py` can drive it the same way.

Conventions
    One TipTorch model holds all N_src = nPointings + nNaturalGS directions: the science pointings first, the NGS
    directions last (as P3 does with getPSDatNGSpositions). PSDs are [N_src, N_wvl, nOtf, nOtf] in nm² on an odd grid
    with the DC at nOtf//2; science PSFs are [nPointings, N_wvl, N_pix, N_pix] normalized to unit sum.
"""

from __future__ import annotations

import ast
import json
from configparser import ConfigParser
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import torch
import yaml
from astropy.io import fits

from mastsel import MavisLO, maskSA, polarToCartesian, psdSetToPsfSet
from tiptorch._config import default_device, default_torch_type
from tiptorch.PSF_models.TipTorch import TipTorch
from tiptorch.managers.config_manager import ConfigManager
from tiptorch.managers.parameter_parser import ParameterParser
from tiptorch.tools.tiptop_integration import (
    RAD_TO_MAS,
    circular_pupil,
    combine_zero_centered_jitters,
    fwhm_to_jitter,
    interpolate_curves,
    pad_PSD_to_even,
    PSF_encircled_energy,
    PSF_ensquared_energy,
    PSF_FWHM,
    PSF_radial_profile,
    tiptilt_covariance_to_jitter,
)

try:
    from ._version import __version__

except ImportError:
    __version__ = 'unknown'

MAX_HEADER_CHARS = 80 # as tiptopUtils.add_hdr_keyword


def _as_list(value, count=None):
    values = np.atleast_1d(value).tolist()
    
    if count is not None and len(values) == 1:
        values = values * count
    if count is not None and len(values) != count:
        raise ValueError(f'Expected {count} values, got {len(values)}')
    
    return values


def _host(value):
    return value.get() if hasattr(value, 'get') else np.asarray(value)


def _add_header_keyword(hdr, section, key, value, i=None, j=None):
    ''' Store one config entry as a HIERARCH card, split over several cards when too long (as tiptopUtils.add_hdr_keyword) '''
    card = ' '.join(str(part) for part in ('HIERARCH ' + section, key, i, j) if part is not None)
    text = str(value)
    
    if len(card) + 4 > MAX_HEADER_CHARS:
        return
    
    while len(card) + len(text) + 5 >= MAX_HEADER_CHARS:
        cut = MAX_HEADER_CHARS - len(card) - 8
        hdr[card + '+'] = text[:cut] + '&&&'
        text = text[cut:]
        
    hdr[card] = text


class baseSimulation:
    """ TipTorch implementation of TIPTOP's `baseSimulation` interface """

    def __init__(self, path, parametersFile, outputDir, outputFile, doConvolve=True, doPlot=False, addSrAndFwhm=True,
                 verbose=False, getHoErrorBreakDown=False, savePSDs=False, ensquaredEnergy=False, eeRadiusInMas=50,
                 backendHO='tiptorch', tiptorch_device=None, tiptorch_dtype=None, tiptorch_kwargs=None):

        self.path, self.parametersFile  = str(path), parametersFile
        self.outputDir, self.outputFile = str(outputDir), outputFile
        self.doConvolve, self.doPlot, self.addSrAndFwhm, self.verbose = doConvolve, doPlot, addSrAndFwhm, verbose
        self.getHoErrorBreakDown, self.savePSDs = getHoErrorBreakDown, savePSDs
        self.ensquaredEnergy, self.eeRadiusInMas = ensquaredEnergy, eeRadiusInMas

        self.firstSimCall = True
        self.doConvolveAsterism = True
        self.pointings_FWHM_mas = None
        self.GFinPSD = False
        self.model = None
        self.results = []

        filename = Path(self.path) / parametersFile
        
        if not filename.suffix:
            filename = filename.with_suffix('.ini')
            if not filename.exists():
                filename = filename.with_suffix('.yml')
                
        self.fullPathFilename = str(filename)
        self.my_data_map = self._read_tiptop_config(filename)
        self._validate_config()

        telescope, science, sensor = self.my_data_map['telescope'], self.my_data_map['sources_science'], self.my_data_map['sensor_science']
        self.tel_radius = float(telescope['TelescopeDiameter']) / 2
        self.wvl = [float(w) for w in _as_list(science['Wavelength'])]
        self.nWvl, self.wvlMax = len(self.wvl), max(self.wvl)
        self.iRef = int(np.argmin(self.wvl)) # P3 evaluates per-source quantities at the shortest wavelength
        self.zenithSrc  = _as_list(science['Zenith'])
        self.azimuthSrc = _as_list(science['Azimuth'], len(self.zenithSrc))
        self.nPointings = len(self.zenithSrc)
        self.pointings  = polarToCartesian(np.array([self.zenithSrc, self.azimuthSrc]))
        self.xxSciencePointigs, self.yySciencePointigs = self.pointings[0, :], self.pointings[1, :]
        self.cartSciencePointingCoords = self.pointings.T.copy()
        self.psInMas = float(sensor['PixelScale'])
        self.nPixPSF = int  (sensor['FieldOfView'])
        self.SupSamp = sensor['Super_Sampling']
        self.jitter_FWHM   = telescope.get('jitter_FWHM')
        self.addFocusError = bool(telescope['glFocusOnNGS'])
        self.LOisOn = 'sensor_LO' in self.my_data_map
        self.nNaturalGS_field = 0
        self.LO_zen_field, self.LO_az_field = [], []

        if self.LOisOn and self.addFocusError and 'sensor_Focus' not in self.my_data_map \
           and max(_as_list(self.my_data_map['sensor_LO']['NumberLenslets'])) == 1:
            raise ValueError('[telescope] glFocusOnNGS is available only if NGS/Focus WFSs have more than one sub-aperture')

        # TipTorch configuration: the LO-independent part is parsed once, the sources are set when the model is built
        self.device = torch.device(tiptorch_device or default_device)
        self.dtype  = tiptorch_dtype or default_torch_type
        self.model_params = self._model_params_from_file(filename)
        
        for key in ('PathPupil', 'PathApodizer', 'PathStaticOn'):
            resolved = self._resolve_data_path(telescope.get(key), filename)
            if resolved is not None:
                self.model_params['telescope'][key] = str(resolved)

        self.tiptorch_kwargs = dict(tiptorch_kwargs or {})
        self.tiptorch_kwargs.setdefault('norm_regime', 'sum')
        self.tiptorch_kwargs.setdefault('oversampling', 1)
        self.tiptorch_kwargs.setdefault('retain_PSDs', getHoErrorBreakDown)
        self.tiptorch_kwargs.setdefault('dtype', self.dtype)
        
        if self.model_params['telescope'].get('PathPupil') is None and 'pupil' not in self.tiptorch_kwargs:
            self.tiptorch_kwargs['pupil'] = circular_pupil(
                int(telescope['Resolution']),
                float(telescope.get('ObscurationRatio', 0)),
                device=self.device,
                dtype=self.dtype
            )
            
        static_path = self.model_params['telescope'].get('PathStaticOn')
        self.static_WFE_nm = None if static_path is None else np.asarray(fits.getdata(static_path), dtype=np.float64)


    # ----------------------------------------- Configuration -----------------------------------------
    def check_section_key(self, primary):
        return primary in self.my_data_map

    def check_config_key(self, primary, secondary):
        return primary in self.my_data_map and secondary in self.my_data_map[primary]


    @staticmethod
    def _read_tiptop_config(filename):
        if not filename.exists():
            raise FileNotFoundError(filename)
        
        if filename.suffix.lower() in ('.yaml', '.yml'):
            with filename.open(encoding='utf-8') as stream:
                return yaml.safe_load(stream)

        parser = ConfigParser(interpolation=None)
        parser.optionxform = str
        parser.read(filename, encoding='utf-8')
        result = {}
        
        for section in parser.sections():
            result[section] = {}
            for key, value in parser.items(section):
                try:
                    result[section][key] = ast.literal_eval(value)
                except (ValueError, SyntaxError):
                    result[section][key] = value.strip()
                    
        return result


    def _validate_config(self):
        ''' The checks and defaults of baseSimulation.__init__ '''
        for section in ('telescope', 'sources_science', 'sensor_science'):
            if section not in self.my_data_map:
                raise ValueError(f"The section '{section}' is missing from the parameter file")
            
        for section, key in (('telescope', 'TelescopeDiameter'), ('sources_science', 'Wavelength')):
            if key not in self.my_data_map[section]:
                raise ValueError(f"'{key}' is missing from section '{section}'")

        science = self.my_data_map['sources_science']
        science.setdefault('Zenith',  [0.0])
        science.setdefault('Azimuth', [0.0])
        
        if len(_as_list(science['Zenith'])) != len(_as_list(science['Azimuth'])):
            raise ValueError("'Zenith' and 'Azimuth' in section 'sources_science' must have the same length")
        
        self.my_data_map['telescope'].setdefault('glFocusOnNGS', False)

        sensor = self.my_data_map['sensor_science']
        super_sampling = sensor.get('Super_Sampling')
        
        if super_sampling is None:
            sensor['Super_Sampling'] = None
        elif isinstance(super_sampling, (int, float)):
            sensor['Super_Sampling'] = [float(super_sampling), 2]  
        elif isinstance(super_sampling, (list, tuple)) and len(super_sampling) in (1, 2):
            sensor['Super_Sampling'] = [float(super_sampling[0]), int(super_sampling[1]) if len(super_sampling) == 2 else 2]
        else:
            raise KeyError('Super_Sampling must be a scalar or list of one/two values.')

        if ('sensor_LO' in self.my_data_map) != ('sources_LO' in self.my_data_map):
            raise KeyError("'sensor_LO' and 'sources_LO' must be defined together")
        
        if 'sensor_LO' in self.my_data_map:
            if 'Wavelength' not in self.my_data_map['sources_LO']:
                raise ValueError("'Wavelength' is missing from section 'sources_LO'")
            
            if 'NumberPhotons' not in self.my_data_map['sensor_LO']:
                raise ValueError("'NumberPhotons' is missing from section 'sensor_LO'")
            
            if 'RTC' not in self.my_data_map or 'SensorFrameRate_LO' not in self.my_data_map['RTC']:
                raise ValueError("'SensorFrameRate_LO' is missing from section 'RTC'")
            
            self.my_data_map['sources_LO'].setdefault('Zenith',  [0.0])
            self.my_data_map['sources_LO'].setdefault('Azimuth', [0.0])


    @staticmethod
    def _resolve_data_path(value, config_file):
        if not isinstance(value, str) or value == '' or value.startswith('$CALIBRATIONS_PATH$'):
            return None # ParameterParser resolves TipTorch calibration tokens itself
        
        path = Path(value)
        candidates = [path] if path.is_absolute() else [base / path for base in (Path.cwd(), *config_file.resolve().parents)]
        for candidate in candidates:
            if candidate.exists():
                return candidate.resolve()
        
        raise FileNotFoundError(f'{value} (from {config_file})')


    def _model_params_from_file(self, filename):
        if filename.suffix.lower() == '.ini':
            return ParameterParser(str(filename)).params
        # TipTorch's parser supplies the defaults its model needs but reads INI only: feed it an INI view of the YAML
        parser = ConfigParser(interpolation=None)
        parser.optionxform = str
        for section, values in self.my_data_map.items():
            if isinstance(values, dict):
                parser[section] = {key: repr(value) for key, value in values.items()}
                
        with TemporaryDirectory() as temp_dir:
            ini_path = Path(temp_dir) / 'tiptop_config.ini'
            with ini_path.open('w', encoding='utf-8') as stream:
                parser.write(stream)
            
            return ParameterParser(str(ini_path)).params


    def configLO(self, astIndex=None):
        source, sensor, RTC = self.my_data_map['sources_LO'], self.my_data_map['sensor_LO'], self.my_data_map['RTC']
        self.LO_zen_field = _as_list(source['Zenith'])
        count = len(self.LO_zen_field)
        self.LO_az_field     = _as_list(source['Azimuth'], count)
        self.LO_wvl          = float(_as_list(source['Wavelength'])[0]) # one wavelength for all the asterism stars
        self.LO_fluxes_field = _as_list(sensor['NumberPhotons'], count)
        self.LO_psInMas      = _as_list(sensor['PixelScale'], count)
        self.LO_freqs_field  = _as_list(RTC['SensorFrameRate_LO'], count)
        self.addLoAlias      = bool(sensor.get('addAliasError', False))

        focus = self.my_data_map.get('sensor_Focus')
        if focus is None:
            self.Focus_fluxes4s_field, self.Focus_psInMas, self.Focus_wvl = self.LO_fluxes_field, self.LO_psInMas, self.LO_wvl
        else:
            self.Focus_fluxes4s_field = _as_list(focus['NumberPhotons'], count)
            self.Focus_psInMas = _as_list(focus['PixelScale'], count)
            self.Focus_wvl = float(_as_list(self.my_data_map.get('sources_Focus', source)['Wavelength'])[0])
            
        self.Focus_freqs_field = _as_list(RTC.get('SensorFrameRate_Focus', self.LO_freqs_field), count)

        self.NGS_fluxes_field   = (np.asarray(self.LO_fluxes_field) * np.asarray(self.LO_freqs_field)).tolist()
        self.Focus_fluxes_field = (np.asarray(self.Focus_fluxes4s_field) * np.asarray(self.Focus_freqs_field)).tolist()
        self.nNaturalGS_field   = count
        self.cartNGSCoords_field = np.array([polarToCartesian(np.array([z, a])) for z, a in zip(self.LO_zen_field, self.LO_az_field)])
        self.currentAsterismIndices = list(range(count))
        self.setAsterismData()


    def setAsterismData(self):
        ids = self.currentAsterismIndices
        self.LO_zen_asterism        = [self.LO_zen_field[i] for i in ids]
        self.LO_az_asterism         = [self.LO_az_field[i] for i in ids]
        self.LO_fluxes_asterism     = [self.LO_fluxes_field[i] for i in ids]
        self.LO_freqs_asterism      = [self.LO_freqs_field[i] for i in ids]
        self.NGS_fluxes_asterism    = [self.NGS_fluxes_field[i] for i in ids]
        self.Focus_fluxes_asterism  = [self.Focus_fluxes_field[i] for i in ids]
        self.cartNGSCoords_asterism = [self.cartNGSCoords_field[i] for i in ids]


    # ----------------------------------------- TipTorch model -----------------------------------------
    def _build_model(self):
        ''' One TipTorch model for the science pointings followed by the NGS directions, sharing one observation '''
        params = deepcopy(self.model_params)
        params['sources_science']['Zenith']  = [float(z) for z in (*self.zenithSrc,  *self.LO_zen_field)]
        params['sources_science']['Azimuth'] = [float(a) for a in (*self.azimuthSrc, *self.LO_az_field)]
        params['NumberSources'] = 1

        manager = ConfigManager()
        manager.select_required_fields(params)
        manager.wrap_scalars_to_lists(params)
        manager.ensure_dimensions(params, 1) # one atmosphere / WFS observation serves all directions
        params['NumberSources'] = len(params['sources_science']['Zenith'])
        self.config_torch = manager.Convert(params, framework='pytorch', device=self.device, dtype=self.dtype)

        model = TipTorch(AO_config=self.config_torch, device=self.device, **self.tiptorch_kwargs)
        self._set_jitter(model, None)

        if self.static_WFE_nm is not None:
            if self.static_WFE_nm.shape != tuple(model.pupil.shape):
                raise ValueError('PathStaticOn map shape does not match the pupil')
            static_WFE = torch.as_tensor(self.static_WFE_nm, device=self.device, dtype=self.dtype)
            phase = 2*torch.pi*1e-9 * static_WFE[None, None] / model.wvl[..., None, None] # [1, N_wvl, N, N]
            model.ComputeStaticOTF(model.pupil * torch.exp(1j*phase))
        return model


    @staticmethod
    def _set_jitter(model, jitter):
        ''' Set the (Jx, Jy, Jxy) jitter of the science pointings; the NGS directions always get zero '''
        zeros = torch.zeros(model.N_src, device=model.device, dtype=model.wvl.dtype)
        if jitter is None:
            model.Jx, model.Jy, model.Jxy = zeros, zeros.clone(), zeros.clone()
        else:
            model.Jx, model.Jy, model.Jxy = (torch.cat([j.to(zeros), zeros[j.numel():]]) for j in jitter)


    def _extra_error_PSD(self):
        ''' P3's extra-error PSDs [N_src, 1, nOtf, nOtf]: science pointings get extraErrorNm, NGS directions extraErrorLoNm (field-dependent) '''
        tel = self.my_data_map['telescope']
        rms_HO, n_NGS = float(tel.get('extraErrorNm', 0)), self.nNaturalGS_field
        HO_shape = (float(tel.get('extraErrorExp', -2)), float(tel.get('extraErrorMin', 0)), float(tel.get('extraErrorMax', 0)))

        RMS_LO = _as_list(tel.get('extraErrorLoNm', rms_HO))
        if len(RMS_LO) == 2: # linear interpolation between the field center and the edge of the technical field
            RMS_LO = np.interp(self.LO_zen_field, [0, float(tel['TechnicalFoV'])/2], RMS_LO).tolist()
        elif len(RMS_LO) == 1:
            RMS_LO = RMS_LO * n_NGS
        else:
            raise ValueError('extraErrorLoNm must be a scalar or [center, edge] values')
        
        LO_shape = (float(tel.get('extraErrorLoExp', HO_shape[0])), float(tel.get('extraErrorLoMin', 0)), float(tel.get('extraErrorLoMax', 0)))

        PSD = 0.0
        if np.sum(RMS_LO) >= 0 and n_NGS > 0:
            PSD = PSD + self.model.ExtraErrorPSD([0.0]*self.nPointings + RMS_LO, *LO_shape)
            rms_HO_per_source = [rms_HO]*self.nPointings + [0.0]*n_NGS
        else:
            rms_HO_per_source = [rms_HO]*(self.nPointings + n_NGS) # negative LO values fall back to the HO extra error
        if rms_HO > 0:
            PSD = PSD + self.model.ExtraErrorPSD(rms_HO_per_source, *HO_shape)
        
        return PSD


    def _prepare_state(self):
        ''' Build the model and compute everything that does not depend on the selected asterism '''
        self.model = self._build_model()
        model = self.model
        
        with torch.no_grad():
            PSD = model.ComputePSD().real.clamp_min(0) # [N_src, N_wvl, nOtf, nOtf] in nm²
            
            if self.LOisOn:
                PSD = PSD * model.TiltFilter() # tip/tilt is handled by the LO loop
                
            elif self.my_data_map['telescope'].get('windPsdFile'):
                wind_path = self._resolve_data_path(self.my_data_map['telescope']['windPsdFile'], Path(self.fullPathFilename))
                PSD = PSD + model.WindShakePSD(fits.getdata(wind_path))
                
            self.PSD_HO = PSD + self._extra_error_PSD()

        if self.getHoErrorBreakDown:
            self.HO_error_budget = model.ErrorBudget(verbose=self.verbose)

        # MASTSEL grid description of the TipTorch PSD, see baseSimulation.doOverallSimulation
        self.N  = int(model.nOtf)
        self.dk = 0.5e9 * float(model.dk) # TIPTOP convention: MASTSEL gets nm² PSDs together with half of the frequency step scaled by 1e9
        self.PSDstep = float(model.dk)
        pixels_per_l_D = self.wvlMax * float(model.rad2mas) / (self.psInMas * 2*self.tel_radius)
        self.overSamp = max(1, int(round(model.sampling_min / pixels_per_l_D)))
        self.pupil = model.pupil.detach().cpu().numpy()

        if self.LOisOn:
            self.ngsPSF()
            self.mLO = MavisLO(self.path, self.parametersFile, verbose=self.verbose)


    # ----------------------------------------- LO sensors (MASTSEL) -----------------------------------------
    def _MASTSEL_PSFs(self, PSD, mask, wavelength, nPixPSF, skip_reshape=False):
        ''' Long-exposure PSFs from TipTorch PSDs [N, nOtf, nOtf] through MASTSEL; odd PSDs are zero-padded to the parity of nPixPSF '''
        PSD = pad_PSD_to_even(PSD) if nPixPSF % 2 == 0 else PSD
        n = PSD.shape[-1]
        freq_range = n * self.PSDstep
        sx = int(2 * np.round(self.tel_radius * freq_range))
        PSFs = psdSetToPsfSet(list(PSD.detach().cpu().numpy()), mask, wavelength, n, sx, 1/self.PSDstep, freq_range, self.dk,
                              nPixPSF, self.wvlMax, self.overSamp, opdMap=self.static_WFE_nm, skip_reshape=skip_reshape)
        images = torch.as_tensor(np.stack([_host(p.sampling) for p in PSFs]), device=self.device, dtype=self.dtype)
        
        return images, float(PSFs[0].pixel_size * RAD_TO_MAS)


    def _sensor_metrics(self, wavelength, pixel_scales, lenslets):
        ''' Strehl, FWHM [mas] and encircled energy of the NGS spots seen by a LO-type sensor, as baseSimulation.ngsPSF '''
        model, n_NGS = self.model, self.nNaturalGS_field
        
        lenslets = _as_list(lenslets)
        if len(lenslets) not in (1, n_NGS):
            raise ValueError('NumberLenslets must be a scalar or have one value per NGS')
        
        mask = maskSA(lenslets, n_NGS, self.pupil)

        with torch.no_grad():
            PSD = self.PSD_HO[self.nPointings:, self.iRef] # [N_NGS, nOtf, nOtf]
            n_sub = torch.as_tensor(lenslets if len(lenslets) == n_NGS else lenslets * n_NGS, device=self.device, dtype=self.dtype).view(-1, 1, 1)
            sub_aperture_filter = model.half_PSD_to_full(model._spatial_filters(model.k, D=model.D/n_sub)[0])
            PSD = torch.where(n_sub > 1, PSD * sub_aperture_filter, PSD)

            psf_scale = self.psInMas * wavelength / self.wvlMax
            skip_reshape = psf_scale / min(pixel_scales) > 1 and self.overSamp > 1
            nPix = self.nPixPSF * self.overSamp if skip_reshape else self.nPixPSF
            images, scale = self._MASTSEL_PSFs(PSD, mask, wavelength, nPix, skip_reshape)

            SR = torch.exp(-PSD.sum(dim=(-2,-1)) * (2*torch.pi*1e-9/wavelength)**2)
            FWHM = torch.sqrt(torch.prod(torch.stack(PSF_FWHM(images, scale)), dim=0))
            radii, EE_curves = PSF_encircled_energy(images, scale)
            EE = interpolate_curves(radii, EE_curves, torch.maximum(FWHM, torch.as_tensor(pixel_scales, device=self.device, dtype=self.dtype)))
            EE = torch.where(2*FWHM >= nPix*scale, torch.ones_like(EE), EE)
            
        return SR.cpu().tolist(), FWHM.cpu().tolist(), EE.cpu().tolist(), mask


    def ngsPSF(self):
        lenslets = self.my_data_map['sensor_LO']['NumberLenslets']
        self.NGS_SR_field, self.NGS_FWHM_mas_field, self.NGS_EE_field, mask = self._sensor_metrics(self.LO_wvl, self.LO_psInMas, lenslets)

        self.NGS_DL_FWHM_mas = None
        
        if self.addLoAlias: # diffraction-limited spot width of every distinct sub-aperture mask
            masks = mask if isinstance(mask, list) else [mask]
            widths = []
            for mask_i in masks:
                zero_PSD = torch.zeros(1, self.N, self.N, device=self.device, dtype=self.dtype)
                images, scale = self._MASTSEL_PSFs(zero_PSD, mask_i, self.LO_wvl, self.nPixPSF)
                widths.append(float(torch.sqrt(torch.prod(torch.stack(PSF_FWHM(images, scale))))))
            self.NGS_DL_FWHM_mas = widths if isinstance(mask, list) else widths[0]

        if self.addFocusError:
            if 'sensor_Focus' in self.my_data_map:
                lenslets = self.my_data_map['sensor_Focus']['NumberLenslets']
                self.Focus_SR_field, self.Focus_FWHM_mas_field, self.Focus_EE_field, _ = self._sensor_metrics(self.Focus_wvl, self.Focus_psInMas, lenslets)
            else:
                self.Focus_SR_field, self.Focus_FWHM_mas_field, self.Focus_EE_field = self.NGS_SR_field, self.NGS_FWHM_mas_field, self.NGS_EE_field


    def _compute_LO(self, astIndex):
        ''' Total tip/tilt (and focus) residual covariances from MavisLO, as the LO part of baseSimulation.doOverallSimulation '''
        field_args = (
            self.cartSciencePointingCoords,
            self.cartNGSCoords_field,
            self.NGS_fluxes_field,
            self.LO_freqs_field,
            self.NGS_SR_field,
            self.NGS_EE_field,
            self.NGS_FWHM_mas_field
        )
        
        focus_args = (
            self.cartNGSCoords_field,
            self.Focus_fluxes_field,
            self.Focus_freqs_field,
            self.Focus_SR_field,
            self.Focus_EE_field,
            self.Focus_FWHM_mas_field
        ) if self.addFocusError else None

        if astIndex is None:
            self.Ctot = _host(self.mLO.computeTotalResidualMatrix(*field_args, aNGS_FWHM_DL_mas=self.NGS_DL_FWHM_mas, doAll=True))
            if self.addFocusError:
                self.CtotFocus = _host(self.mLO.computeFocusTotalResidualMatrix(*focus_args))
                self.GF_res = float(np.sqrt(max(self.CtotFocus[0], 0)))
                self.PSD = self.PSD + self.model.FocusErrorPSD(self.GF_res)
                self.GFinPSD = True
        else:
            if self.firstSimCall:
                self.mLO.computeTotalResidualMatrix(*field_args, aNGS_FWHM_DL_mas=self.NGS_DL_FWHM_mas, doAll=False)
                
                if self.addFocusError:
                    self.mLO.computeFocusTotalResidualMatrix(*focus_args)
            
            # Discard guide stars with less than 1 photon per frame per sub-aperture
            if np.min(self.NGS_fluxes_asterism) < 1 and np.max(self.NGS_fluxes_asterism) > 1:
                keep = [i for i, flux in enumerate(self.NGS_fluxes_asterism) if flux > 1]
                self.currentAsterismIndices = [self.currentAsterismIndices[i] for i in keep]
                self.setAsterismData()
                
            self.Ctot = _host(self.mLO.computeTotalResidualMatrixI(self.currentAsterismIndices, self.cartSciencePointingCoords,
                                                                   np.asarray(self.cartNGSCoords_asterism), self.NGS_fluxes_asterism))
            if self.addFocusError:
                self.CtotFocus = _host(self.mLO.computeFocusTotalResidualMatrixI(self.currentAsterismIndices, np.asarray(self.cartNGSCoords_asterism),
                                                                                 self.Focus_fluxes_asterism))
                self.GF_res = float(np.sqrt(max(self.CtotFocus[0], 0)))
                self.GFinPSD = False

        self.LO_res = np.sqrt(np.maximum(np.trace(self.Ctot, axis1=-2, axis2=-1), 0))
        covariance = torch.as_tensor(self.Ctot, device=self.device, dtype=self.dtype)
        self.LO_jitter = tiptilt_covariance_to_jitter(covariance, 2*self.tel_radius) # (Jx, Jy [mas], Jxy [deg]) per pointing
        Jx, Jy, Jxy = (j.detach().cpu().numpy() for j in self.LO_jitter)
        self.cov_ellipses = np.column_stack((np.deg2rad(Jxy), Jx, Jy)) # (angle [rad], sigma_max, sigma_min [mas]) as MavisLO.ellipsesFromCovMats


    # ----------------------------------------- Science PSFs -----------------------------------------
    def _combined_jitter(self):
        ''' Science jitter kernel: the LO residual convolved with the optional telescope jitter, as baseSimulation.finalConvolution '''
        telescope = fwhm_to_jitter(self.jitter_FWHM, device=self.device, dtype=self.dtype, count=self.nPointings) if self.jitter_FWHM is not None else None
        LO = self.LO_jitter if self.LOisOn else None
        if LO is None:
            return telescope
        if telescope is None:
            return LO
        return combine_zero_centered_jitters(*LO, *telescope)


    def _pixel_centered(self, PSFs):
        ''' TipTorch pads its odd OTF grid to an even FieldOfView on the right/bottom, which leaves the PSF peak on pixel N//2 - 1; P3 and MASTSEL center even PSFs on pixel N//2 '''
        if self.nPixPSF % 2 == 0 and self.model.nOtf < self.nPixPSF:
            return torch.roll(PSFs, (1, 1), dims=(-2, -1))
        return PSFs


    def _science_PSFs(self, jitter):
        ''' Science PSFs [nPointings, N_wvl, N_pix, N_pix] of the current PSD with the given jitter '''
        self._set_jitter(self.model, jitter)
        with torch.no_grad():
            return self._pixel_centered(self.model(PSD=self.PSD)[:self.nPointings])


    def _FWHM_per_wavelength(self, PSFs):
        ''' Mean FWHM [mas] of PSFs [nPointings, N_wvl, N, N] as nested lists [N_wvl][nPointings], or [nPointings] for one wavelength '''
        FWHM = 0.5 * sum(PSF_FWHM(PSFs, self.psInMas)) # as P3's getFWHM(nargout=1)
        FWHM = FWHM.T.cpu().tolist()
        return FWHM if self.nWvl > 1 else FWHM[0]


    def finalPSF(self, astIndex=None):
        if astIndex is None or self.firstSimCall:
            self.pointings_FWHM_mas = self._FWHM_per_wavelength(self._science_PSFs(None)) # jitter-free HO PSFs

        if self.LOisOn and not self.doConvolveAsterism:
            self.results = [] # asterism selection only needs the covariance ellipses
            return

        PSFs = self._science_PSFs(self._combined_jitter() if (not self.LOisOn or self.doConvolve) else None)
        self.PSF_tensor = PSFs
        self.results = PSFs
        PSFs = PSFs.cpu().numpy()
        self.cubeResultsArray = np.moveaxis(PSFs, 1, 0) if self.nWvl > 1 else PSFs[:, 0]
        self.cubeResults = self.cubeResultsArray


    def doOverallSimulation(self, astIndex=None):
        if self.LOisOn:
            self.configLO(astIndex)
        self.results = []

        if astIndex is None or self.firstSimCall:
            self._prepare_state()

        self.PSD = self.PSD_HO
        if self.LOisOn:
            self._compute_LO(astIndex)
        self.HO_res = self.PSD[:self.nPointings, self.iRef].sum(dim=(-2,-1)).clamp_min(0).sqrt().cpu().numpy()
        self.finalPSF(astIndex)
        self.firstSimCall = False

        if astIndex is None:
            with torch.no_grad():
                self.PSF_DL = self._pixel_centered(self.model.DLPSF()[0])                                          # [N_wvl, N_pix, N_pix]
                self.PSF_OL = self._pixel_centered(self.model.OLPSF(include_static=True, include_jitter=False)[0]) # [N_wvl, N_pix, N_pix]
            self.psfDL = self.PSF_DL[self.iRef].cpu().numpy()
            self.psfOL = self.PSF_OL[self.iRef].cpu().numpy()
            self.computePSF1D()
            if self.verbose:
                print('HO_res [nm]:', self.HO_res)
                if self.LOisOn:
                    print('LO_res [nm]:', self.LO_res)
                if hasattr(self, 'GF_res'):
                    print('GF_res [nm]:', self.GF_res)

        if self.doPlot and len(self.results):
            import matplotlib.pyplot as plt
            plt.imshow(np.log10(np.abs(self.cubeResultsArray[0] if self.nWvl == 1 else self.cubeResultsArray[0, 0]) + 1e-12))
            plt.show()


    # ----------------------------------------- Metrics and outputs -----------------------------------------
    def computeMetrics(self):
        HO_res = np.asarray(self.HO_res, dtype=float)
        if not len(self.results): # covariance-only asterism evaluation, as baseSimulation.computeMetrics
            variance = HO_res**2 + (self.LO_res**2 if self.LOisOn else 0) + (self.GF_res**2 if self.addFocusError and not self.GFinPSD else 0)
            self.penalty = [float(np.sqrt(np.mean(variance)))]
            self.sr = [float(np.exp(-4*np.pi**2 * self.penalty[0]**2 / (self.wvl[self.iRef]*1e9)**2))]
            scale = np.pi/(180*3600*1000) * 2*self.tel_radius / 4e-9
            FWHM_LO = 2.355 * self.LO_res / scale / np.sqrt(2) if self.LOisOn else 0.0
            self.fwhm = [np.sqrt(FWHM_LO**2 + np.asarray(self.pointings_FWHM_mas)**2)]
            self.ee = [0]
            return

        variance = HO_res**2 + (self.LO_res**2 if self.LOisOn else 0) + (self.GF_res**2 if self.LOisOn and self.addFocusError and not self.GFinPSD else 0)
        self.penalty = np.sqrt(variance).tolist()

        with torch.no_grad():
            PSFs = self.PSF_tensor # [nPointings, N_wvl, N, N]
            peak = lambda x: x.amax(dim=(-2,-1)) / x.sum(dim=(-2,-1))
            SR = peak(PSFs) / peak(self.PSF_DL)
            FWHM = 0.5 * sum(PSF_FWHM(PSFs, self.psInMas))
            if self.ensquaredEnergy:
                EE_curves = PSF_ensquared_energy(PSFs)
                radii = (torch.arange(EE_curves.shape[-1], device=self.device, dtype=self.dtype) + 0.5) * self.psInMas
            else:
                radii, EE_curves = PSF_encircled_energy(PSFs, self.psInMas)
            EE = interpolate_curves(radii, EE_curves, self.eeRadiusInMas)

        to_lists = lambda x: x.T.cpu().tolist() if self.nWvl > 1 else x[:, 0].cpu().tolist() # [N_wvl][nPointings] or [nPointings]
        self.sr, self.fwhm, self.ee = to_lists(SR), to_lists(FWHM), to_lists(EE)


    def computePSF1D(self):
        ''' Radial profiles about the PSF peaks, as baseSimulation.computePSF1D without the Super_Sampling interpolation '''
        with torch.no_grad():
            radii, profiles = PSF_radial_profile(self.PSF_tensor, self.psInMas) # [nPointings, N_wvl, n_bins]
            keep = radii < self.psInMas * self.nPixPSF / 2
            self.psf1d_radius = radii[keep].cpu().numpy()
            self.psf1d = profiles[..., keep].permute(1, 0, 2).cpu().numpy() # [N_wvl, nPointings, n_bins]
        self.psf1d_radius_list_list = np.broadcast_to(self.psf1d_radius, self.psf1d.shape).copy()
        self.psf1d_data = np.vstack((self.psf1d_radius_list_list, self.psf1d))


    def savePSFprofileJSON(self):
        output = Path(self.outputDir) / f'{self.outputFile}1D_PSF.json'
        output.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            'execution_infos': {'TIME': datetime.now().strftime('%Y%m%d_%H%M%S'), 'TIPTOP version': __version__},
            'infos': self.my_data_map,
            'psf': {'radius': self.psf1d_radius.tolist(), 'psf': self.psf1d.tolist()},
        }
        output.write_text(json.dumps(payload, default=str), encoding='utf-8')


    def saveResults(self):
        ''' FITS output with the HDUs and header keywords of baseSimulation.saveResults '''
        if not hasattr(self, 'psfOL'):
            raise RuntimeError('Run doOverallSimulation() without astIndex before saving')
        
        now = datetime.now().strftime('%Y%m%d_%H%M%S')

        hdus = [fits.PrimaryHDU(), fits.ImageHDU(data=self.cubeResultsArray), fits.ImageHDU(data=self.psfOL), fits.ImageHDU(data=self.psfDL)]
        
        if self.savePSDs:
            hdus.append(fits.ImageHDU(data=self.PSD[:self.nPointings, self.iRef].cpu().numpy())) # [nPointings, nOtf, nOtf] as baseSimulation
       
        hdus.append(fits.ImageHDU(data=self.psf1d_data))

        hdr0 = hdus[0].header
        hdr0['TIME'], hdr0['TIPTOP_V'] = now, __version__
        for section, entries in self.my_data_map.items():
            for key, value in entries.items():
                if isinstance(value, list):
                    for i, item in enumerate(value):
                        if isinstance(item, list):
                            for j, element in enumerate(item):
                                _add_header_keyword(hdr0, section, key, element, i, j)
                        else:
                            _add_header_keyword(hdr0, section, key, item, i)
                else:
                    _add_header_keyword(hdr0, section, key, value)

        hdr1 = hdus[1].header
        hdr1['TIME'], hdr1['CONTENT'], hdr1['SIZE'] = now, 'PSF CUBE', str(self.cubeResultsArray.shape)
        
        if self.nWvl > 1:
            for i, wavelength in enumerate(self.wvl):
                hdr1['WL_NM' + str(i).zfill(3)] = str(int(wavelength*1e9))
        else:
            hdr1['WL_NM'] = str(int(self.wvl[0]*1e9))
       
        hdr1['PIX_MAS'] = str(self.psInMas)
        hdr1['CC'] = f'CARTESIAN COORD. IN ASEC OF THE {self.nPointings} SOURCES'
        
        for i in range(self.nPointings):
            hdr1['CCX' + str(i).zfill(4)] = round(float(self.pointings[0, i]), 3)
            hdr1['CCY' + str(i).zfill(4)] = round(float(self.pointings[1, i]), 3)
        
        hdr1['RESH'] = 'High Order residual in nm RMS'
        
        for i, residual in enumerate(self.HO_res):
            hdr1['RESH' + str(i).zfill(4)] = round(float(residual), 3)
        
        if self.LOisOn:
            hdr1['RESL'] = 'Low Order residual in nm RMS'
            for i, residual in enumerate(self.LO_res):
                hdr1['RESL' + str(i).zfill(4)] = round(float(residual), 3)
        if hasattr(self, 'GF_res'):
            hdr1['RESF'] = 'Global Focus residual in nm RMS (included in PSD)'
            hdr1['RESF0000'] = round(float(self.GF_res), 3)

        if self.addSrAndFwhm:
            self.computeMetrics()
            for i in range(self.nWvl):
                multi = self.nWvl > 1
                sr, fwhm, ee = (metric[i] if multi else metric for metric in (self.sr, self.fwhm, self.ee))
                wTxt, fTxt, eTxt, Nfill = ('W' + str(i).zfill(2), 'FW', 'EE', 2) if multi else ('', 'FWHM', 'EE' + str(np.round(self.eeRadiusInMas)), 4)
                for j in range(self.nPointings):
                    hdr1['SR' + str(j).zfill(Nfill) + wTxt] = round(float(sr[j]), 5)
                for j in range(self.nPointings):
                    hdr1[fTxt + str(j).zfill(Nfill) + wTxt] = round(float(fwhm[j]), 3)
                for j in range(self.nPointings):
                    hdr1[eTxt + str(j).zfill(Nfill) + wTxt] = round(float(ee[j]), 5)

        for hdu, content in zip(hdus[2:], ('OPEN-LOOP PSF', 'DIFFRACTION LIMITED PSF', *(('High Order PSD',) if self.savePSDs else ()), 'Final PSFs profiles')):
            hdu.header['TIME'], hdu.header['CONTENT'], hdu.header['SIZE'] = now, content, str(hdu.data.shape)
        
        if self.SupSamp:
            hdus[-1].header['SAMP_MAS'] = str(self.SupSamp[0])

        output = Path(self.outputDir) / f'{self.outputFile}.fits'
        output.parent.mkdir(parents=True, exist_ok=True)
        fits.HDUList(hdus).writeto(output, overwrite=True)
        
        if self.verbose:
            print('Output cube shape:', self.cubeResultsArray.shape)
