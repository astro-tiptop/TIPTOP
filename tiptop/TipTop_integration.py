"""TIPTOP simulation backed by TipTorch for high-order science PSFs.

MASTSEL still supplies the low-order covariance and guide-star sensor PSFs.
This module is deliberately independent of P3 so it can be imported directly
while backend selection in TIPTOP is developed separately.
"""

from __future__ import annotations

import ast
from configparser import ConfigParser
from copy import deepcopy
from datetime import datetime
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from astropy.io import fits
import numpy as np
from scipy.special import j1
import torch
import yaml

from mastsel import MavisLO, maskSA, polarToCartesian, psdSetToPsfSet
from tiptorch._config import default_device, default_torch_type
from tiptorch.PSF_models.TipTorch import TipTorch
from tiptorch.managers.config_manager import ConfigManager
from tiptorch.managers.parameter_parser import ParameterParser
from tiptorch.tools.tiptop_integration import (
    RAD_TO_MAS, circular_pupil, combine_zero_centered_jitters,
    extra_error_psd, fwhm_to_jitter,
    psf_encircled_energy, psf_ensquared_energy, psf_fwhm,
    tiptilt_covariance_to_jitter, tiptilt_rejection_filter,
    wind_shake_psd,
)


def _as_list(value, count=None):
    values = np.atleast_1d(value).tolist()
    
    if count is not None and len(values) == 1:
        values *= count
        
    if count is not None and len(values) != count:
        raise ValueError(f"Expected {count} values, got {len(values)}")
    
    return values


def _host(value):
    return value.get() if hasattr(value, 'get') else np.asarray(value)


def _piston_filter(diameter, size, dk):
    """Circular subaperture piston filter on a centered PSD grid."""
    frequencies = (np.arange(size) - size // 2) * dk
    fx, fy = np.meshgrid(frequencies, frequencies)
    x = np.pi * diameter * np.hypot(fx, fy)
    ratio = np.full_like(x, 0.5)
    np.divide(j1(x), x, out=ratio, where=x != 0)
    return np.maximum(1.0 - 4.0 * ratio**2, 0.0)


class baseSimulation:
    """TipTorch implementation of TIPTOP's ``baseSimulation`` interface."""

    def __init__(self, path, parametersFile, outputDir, outputFile,
                 doConvolve=True, doPlot=False, addSrAndFwhm=True,
                 verbose=False, getHoErrorBreakDown=False, savePSDs=False,
                 ensquaredEnergy=False, eeRadiusInMas=50,
                 tiptorch_device=None, tiptorch_dtype=None, tiptorch_kwargs=None):
        
        self.path = str(path)
        self.parametersFile = parametersFile
        self.outputDir = str(outputDir)
        self.outputFile = outputFile
        self.doConvolve = doConvolve
        self.doPlot = doPlot
        self.addSrAndFwhm = addSrAndFwhm
        self.verbose = verbose
        self.getHoErrorBreakDown = getHoErrorBreakDown
        self.savePSDs = savePSDs
        self.ensquaredEnergy = ensquaredEnergy
        self.eeRadiusInMas = eeRadiusInMas
        self.firstSimCall = True
        self.doConvolveAsterism = True
        self.pointings_FWHM_mas = None
        self.GFinPSD = False
        self.LO_tiptilt_cov = None

        filename = Path(self.path) / parametersFile
        if not filename.suffix:
            filename = filename.with_suffix('.ini')
            
            if not filename.exists():
                filename = filename.with_suffix('.yml')
                
        self.fullPathFilename = str(filename)
        self.my_data_map = self._read_tiptop_config(filename)
        self._validate_config()

        manager = ConfigManager()
        source_count = len(_as_list(self.my_data_map['sources_science'].get('Zenith', [0.0])))
        model_config = self._model_config_from_file(filename)
        model_config['NumberSources'] = 1
        manager.select_required_fields(model_config)
        manager.wrap_scalars_to_lists(model_config)
        # One atmospheric/WFS observation can serve several science directions.
        manager.ensure_dimensions(model_config, 1)
        
        model_config['NumberSources'] = source_count
        for key in ('PathPupil', 'PathApodizer', 'PathStaticOn'):
            raw_path = self.my_data_map['telescope'].get(key)
            resolved = self._resolve_data_path(raw_path, filename)
            if resolved is not None:
                model_config['telescope'][key] = str(resolved)
                
        device = torch.device(tiptorch_device or default_device)
        dtype = tiptorch_dtype or default_torch_type
        self.config_torch = manager.Convert( model_config, framework='pytorch', device=device, dtype=dtype, )
        
        kwargs = dict(tiptorch_kwargs or {})
        kwargs.setdefault('norm_regime', 'sum')
        kwargs.setdefault('oversampling', 1)
        kwargs.setdefault('retain_PSDs', True)
        kwargs.setdefault('dtype', dtype)
        
        if model_config['telescope'].get('PathPupil') is None and 'pupil' not in kwargs:
            telescope = self.my_data_map['telescope']
            kwargs['pupil'] = circular_pupil(
                int(telescope['Resolution']),
                float(telescope.get('ObscurationRatio', 0)),
                device=device, dtype=dtype)
            
        self.model = TipTorch(AO_config=self.config_torch, device=device, **kwargs)
        self._static_wfe = None
        static_path = model_config['telescope'].get('PathStaticOn')
        
        if static_path is not None:
            static_map = np.asarray(fits.getdata(static_path), dtype=np.float64)
            if static_map.shape != tuple(self.model.pupil.shape):
                raise ValueError('PathStaticOn map shape does not match the pupil')
            
            if not np.isfinite(static_map).all():
                raise ValueError('PathStaticOn map contains non-finite values')
            
            self._static_wfe = torch.as_tensor(static_map, device=device, dtype=dtype)
            phase = (2 * torch.pi * 1e-9 * self._static_wfe[None, None] / self.model.wvl[..., None, None])
            self.model.ComputeStaticOTF(self.model.pupil * torch.exp(1j * phase))
            
        self.nPointings = int(self.model.N_src)
        self.nPixPSF = int(self.model.N_pix)
        self.psInMas = float(torch.as_tensor(self.model.psInMas).item())
        self.overSamp = int(self.model.oversampling)
        self.wvl = self.model.wvl.detach().cpu().flatten().tolist()
        self.wvlMax = max(self.wvl)
        self.nWvl = len(self.wvl)
        self.tel_radius = float(self.my_data_map['telescope']['TelescopeDiameter']) / 2
        science = self.my_data_map['sources_science']
        self.zenithSrc = _as_list(science.get('Zenith', [0.0]))
        self.azimuthSrc = _as_list(science.get('Azimuth', [0.0]), len(self.zenithSrc))
        self.pointings = polarToCartesian(np.array([self.zenithSrc, self.azimuthSrc]))
        self.cartSciencePointingCoords = self.pointings.T.copy()
        self.SupSamp = self.my_data_map['sensor_science'].get('Super_Sampling')
        self.jitter_FWHM = self.my_data_map['telescope'].get('jitter_FWHM')
        self.addFocusError = bool(self.my_data_map['telescope'].get('glFocusOnNGS', False))
        self.LOisOn = 'sensor_LO' in self.my_data_map
        self.nNaturalGS_field = 0
        
        if self.LOisOn and self.addFocusError and 'sensor_Focus' not in self.my_data_map:
            lenslets = _as_list(self.my_data_map['sensor_LO'].get('NumberLenslets', [1]))
            if max(lenslets) == 1:
                raise ValueError('NGS focus correction needs more than one LO lenslet')


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


    @staticmethod
    def _resolve_data_path(value, config_file):
        if value in (None, '') or not isinstance(value, str):
            return None
        path = Path(value)
        if path.is_absolute():
            if not path.exists():
                raise FileNotFoundError(path)
            return path
        for base in (Path.cwd(), *config_file.resolve().parents):
            candidate = base / path
            if candidate.exists():
                return candidate.resolve()
        if value.startswith('$CALIBRATIONS_PATH$'):
            return None  # ParameterParser may already have resolved this token.
        raise FileNotFoundError(f'{value} (from {config_file})')


    def _model_config_from_file(self, filename):
        if filename.suffix.lower() == '.ini':
            return ParameterParser(str(filename)).params
        
        # TipTorch's parser supplies the defaults required by its model but
        # currently reads INI only. Feed it a temporary INI view of the YAML.
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


    def _validate_config(self):
        for section in ('telescope', 'sources_science', 'sensor_science'):
            if section not in self.my_data_map:
                raise ValueError(f"Missing [{section}] section")
            
        if ('sensor_LO' in self.my_data_map) != ('sources_LO' in self.my_data_map):
            raise ValueError('sensor_LO and sources_LO must be defined together')
        
        if 'sensor_LO' in self.my_data_map:
            if 'RTC' not in self.my_data_map or 'SensorFrameRate_LO' not in self.my_data_map['RTC']:
                raise ValueError('RTC.SensorFrameRate_LO is required for LO sensing')


    def configLO(self, astIndex=None):
        source = self.my_data_map['sources_LO']
        sensor = self.my_data_map['sensor_LO']
        rtc = self.my_data_map['RTC']
        self.LO_zen_field = _as_list(source.get('Zenith', [0.0]))
        
        count = len(self.LO_zen_field)
        self.LO_az_field     = _as_list(source.get('Azimuth', [0.0]), count)
        self.LO_wvl          = float(_as_list(source['Wavelength'])[0])
        self.LO_fluxes_field = _as_list(sensor['NumberPhotons'], count)
        self.LO_psInMas      = _as_list(sensor['PixelScale'], count)
        self.LO_freqs_field  = _as_list(rtc['SensorFrameRate_LO'], count)
        
        self.addLOAlias = bool(sensor.get('addAliasError', False))
        focus = self.my_data_map.get('sensor_Focus')
        
        if focus is None:
            self.Focus_fluxes4s_field = self.LO_fluxes_field
            self.Focus_psInMas = self.LO_psInMas
            self.Focus_wvl = self.LO_wvl
        else:
            self.Focus_fluxes4s_field = _as_list(focus['NumberPhotons'], count)
            self.Focus_psInMas = _as_list(focus['PixelScale'], count)
            focus_source = self.my_data_map.get('sources_Focus', source)
            self.Focus_wvl = float(_as_list(focus_source['Wavelength'])[0])
            
        self.Focus_freqs_field = _as_list(rtc.get('SensorFrameRate_Focus', self.LO_freqs_field), count)
        self.NGS_fluxes_field = (np.asarray(self.LO_fluxes_field) * np.asarray(self.LO_freqs_field)).tolist()
        self.Focus_fluxes_field = (np.asarray(self.Focus_fluxes4s_field) * np.asarray(self.Focus_freqs_field)).tolist()
        self.nNaturalGS_field = count
        
        self.cartNGSCoords_field = np.array([ polarToCartesian(np.array([z, a])) for z, a in zip(self.LO_zen_field, self.LO_az_field) ])
        if not hasattr(self, 'currentAsterismIndices') or astIndex is None:
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
        self.cartNGSCoords_asterism = self.cartNGSCoords_field[ids].tolist()

    def _source_model(self, zenith, azimuth, wavelength):
        config = deepcopy(self.config_torch)
        dtype, device = self.model.wvl.dtype, self.model.device
        config['NumberSources'] = len(zenith)
        config['sources_science']['Zenith']     = torch.as_tensor(zenith, device=device, dtype=dtype)
        config['sources_science']['Azimuth']    = torch.as_tensor(azimuth, device=device, dtype=dtype)
        config['sources_science']['Wavelength'] = torch.as_tensor([wavelength], device=device, dtype=dtype)
        
        model = TipTorch(
            AO_config=config,
            device=device, 
            pupil=self.model.pupil,
            dtype=self.model.dtype,
            norm_regime='sum',
            oversampling=self.overSamp
        )
        
        self._zero_jitter(model)
        return model


    def _extra_error(self, model, *, LO=False):
        telescope = self.my_data_map['telescope']
        if LO:
            amplitudes = _as_list(telescope.get(
                'extraErrorLoNm', telescope.get('extraErrorNm', 0)))
            if len(amplitudes) == 2:
                field_radius = float(telescope['TechnicalFoV']) / 2
                amplitudes = np.interp(
                    self.LO_zen_field, [0, field_radius], amplitudes).tolist()
            elif len(amplitudes) == 1:
                amplitudes *= self.nNaturalGS_field
            else:
                raise ValueError('extraErrorLoNm must be a scalar or two field-edge values')
            if all(value < 0 for value in amplitudes):
                amplitudes = [float(telescope.get('extraErrorNm', 0))] * len(amplitudes)
            
            exponent = float(telescope.get('extraErrorLoExp', telescope.get('extraErrorExp', -2)))
            minimum  = float(telescope.get('extraErrorLoMin', telescope.get('extraErrorMin', 0)))
            maximum  = float(telescope.get('extraErrorLoMax', telescope.get('extraErrorMax', 0)))
            
        else:
            amplitudes = [float(telescope.get('extraErrorNm', 0))]
            exponent   = float(telescope.get('extraErrorExp', -2))
            minimum    = float(telescope.get('extraErrorMin', 0))
            maximum    = float(telescope.get('extraErrorMax', 0))
            
        if not any(value > 0 for value in amplitudes):
            return None
        
        shape = extra_error_psd(model, 1.0, exponent, minimum, maximum)
        values = torch.as_tensor(amplitudes, device=model.device, dtype=model.wvl.dtype)
        return values[:, None, None, None].square() * shape


    @staticmethod
    def _zero_jitter(model):
        model.Jx  = torch.zeros(model.N_src, device=model.device, dtype=model.wvl.dtype)
        model.Jy  = torch.zeros_like(model.Jx)
        model.Jxy = torch.zeros_like(model.Jx)


    def _sensor_PSFs(self, wavelength, pixel_scales, lenslets):
        """Use MASTSEL to make subaperture PSFs from TipTorch's NGS PSDs."""
        sensor_model = self._source_model(
            self.LO_zen_field, self.LO_az_field, wavelength)
    
        with torch.no_grad():
            psd_tensor = sensor_model.ComputePSD().detach().real.clamp_min(0)
            psd_tensor = psd_tensor * tiptilt_rejection_filter(sensor_model)
            extra = self._extra_error(sensor_model, LO=True)
            if extra is not None:
                psd_tensor = psd_tensor + extra
                
            psd = psd_tensor.cpu().numpy()[:, 0].copy()
            
        n = psd.shape[-1]
        dk = float(sensor_model.dk)
        n_lenslets = _as_list(lenslets)
        
        if len(n_lenslets) not in (1, self.nNaturalGS_field):
            raise ValueError('NumberLenslets must be scalar or one value per NGS')
        
        for i in range(self.nNaturalGS_field):
            lenslets_i = n_lenslets[i] if len(n_lenslets) > 1 else n_lenslets[0]
            if lenslets_i != 1:
                psd[i] *= _piston_filter(2 * self.tel_radius / lenslets_i, n, dk)
                
        pupil = sensor_model.pupil.detach().cpu().numpy()
        mask = maskSA(n_lenslets, self.nNaturalGS_field, pupil)
        psf_scale = self.psInMas * wavelength / self.wvlMax
        skip_reshape = psf_scale > min(pixel_scales) and self.overSamp > 1
        size = self.nPixPSF * self.overSamp if skip_reshape else self.nPixPSF
        
        if skip_reshape:
            psf_scale /= self.overSamp
            
        if size % 2 != n % 2 and (self.overSamp == 1 or skip_reshape):
            size += 1
            
        psfs = psdSetToPsfSet(
            psd, mask, wavelength, n, int(round(2 * self.tel_radius * n * dk)),
            1 / dk, n * dk, 1e9 * dk, size, self.wvlMax, self.overSamp,
            opdMap=(self._static_wfe.cpu().numpy() if self._static_wfe is not None else None), skip_reshape=skip_reshape)
        
        return psd, psfs, mask, float(psfs[0].pixel_size * RAD_TO_MAS)


    def _sensor_metrics(self, wavelength, pixel_scales, lenslets):
        psd, psfs, mask, scale = self._sensor_PSFs(wavelength, pixel_scales, lenslets)
        sr, fwhm, ee = [], [], []
        
        for i, image in enumerate(psfs):
            array = _host(image.sampling)
            sr.append(float(np.exp(-psd[i].sum() * (2 * np.pi * 1e-9 / wavelength)**2)))
            fx, fy = psf_fwhm(array, scale)
            width = float(np.sqrt(fx * fy))
            fwhm.append(width)
            
            if 2 * width >= array.shape[-1] * scale:
                ee.append(1.0)
            else:
                curve, radii = psf_encircled_energy(array, scale)
                ee.append(float(np.interp(max(width, pixel_scales[i]), radii, curve)))
        return sr, fwhm, ee, mask, psd, scale


    def ngsPSF(self):
        lenslets = self.my_data_map['sensor_LO'].get('NumberLenslets', [1])
        (self.NGS_SR_field, self.NGS_FWHM_mas_field, self.NGS_EE_field,mask, psd, scale) = self._sensor_metrics(self.LO_wvl, self.LO_psInMas, lenslets)
         
        self.NGS_DL_FWHM_mas = None
        if self.addLOAlias:
            masks = mask if isinstance(mask, list) else [mask]
            widths = []
            n = psd.shape[-1]
            dk = float(self._source_model(
                self.LO_zen_field, self.LO_az_field, self.LO_wvl).dk)
            for pupil in masks:
                size = self.nPixPSF + (self.nPixPSF % 2 != n % 2)
                # psdSetToPsfSet handles the same grid and sizing as the NGS PSFs.
                dl = psdSetToPsfSet(
                    [np.zeros_like(psd[0])], pupil, self.LO_wvl, n,
                    int(round(2 * self.tel_radius * n * dk)), 1 / dk,
                    n * dk, 1e9 * dk, size, self.wvlMax,
                    self.overSamp)[0]
                fx, fy = psf_fwhm(_host(dl.sampling), dl.pixel_size * RAD_TO_MAS)
                widths.append(float(np.sqrt(fx * fy)))
            self.NGS_DL_FWHM_mas = widths if isinstance(mask, list) else widths[0]
            
        if not self.addFocusError:
            return
        
        if 'sensor_Focus' not in self.my_data_map:
            self.Focus_SR_field = self.NGS_SR_field
            self.Focus_FWHM_mas_field = self.NGS_FWHM_mas_field
            self.Focus_EE_field = self.NGS_EE_field
            
        else:
            lenslets = self.my_data_map['sensor_Focus'].get('NumberLenslets', [1])
            (self.Focus_SR_field, self.Focus_FWHM_mas_field, self.Focus_EE_field, _, _, _) = self._sensor_metrics(self.Focus_wvl, self.Focus_psInMas, lenslets)


    def _focus_filter(self):
        half = self.model.FocusFilter(self.model.k2).real
        full = self.model.half_PSD_to_full(half).real
        return full / full.sum(dim=(-2, -1), keepdim=True).clamp_min(1e-30)


    def _prepare_state(self):
        self._zero_jitter(self.model)
        with torch.no_grad():
            self.PSD = self.model.ComputePSD().detach().real.clamp_min(0).clone()
            if self.LOisOn:
                self.PSD *= tiptilt_rejection_filter(self.model)
            extra = self._extra_error(self.model)
            if extra is not None:
                self.PSD += extra
            wind_file = self.my_data_map['telescope'].get('windPsdFile')
            if not self.LOisOn and wind_file:
                wind_path = self._resolve_data_path(
                    wind_file, Path(self.fullPathFilename))
                rtc = self.my_data_map['RTC']
                self.PSD += wind_shake_psd(
                    self.model, fits.getdata(wind_path),
                    float(_as_list(rtc['SensorFrameRate_HO'])[0]),
                    float(_as_list(rtc['LoopGain_HO'])[0]),
                    float(_as_list(rtc['LoopDelaySteps_HO'])[0]))
        self._PSD_HO = self.PSD.clone()
        self._ho_fwhm = None
        self.HO_res = self.PSD.sum(dim=(-2, -1)).clamp_min(0).sqrt().cpu().numpy()[:, 0]
        
        if self.LOisOn:
            self.ngsPSF()
            self.mLO = MavisLO(self.path, self.parametersFile, verbose=self.verbose)


    def _compute_LO(self, astIndex):
        if not self.LOisOn:
            return
        
        args = (self.cartSciencePointingCoords, self.cartNGSCoords_field,
                self.NGS_fluxes_field, self.LO_freqs_field,
                self.NGS_SR_field, self.NGS_EE_field, self.NGS_FWHM_mas_field)
        
        if astIndex is None:
            self.LO_tiptilt_cov = _host(self.mLO.computeTotalResidualMatrix(*args, aNGS_FWHM_DL_mas=self.NGS_DL_FWHM_mas, doAll=True))
            
        else:
            if self.firstSimCall:
                self.mLO.computeTotalResidualMatrix(*args, aNGS_FWHM_DL_mas=self.NGS_DL_FWHM_mas, doAll=False)
                
                if self.addFocusError:
                    self.mLO.computeFocusTotalResidualMatrix(
                        self.cartNGSCoords_field, self.Focus_fluxes_field,
                        self.Focus_freqs_field, self.Focus_SR_field,
                        self.Focus_EE_field, self.Focus_FWHM_mas_field)
                    
            keep = [j for j, f in enumerate(self.NGS_fluxes_asterism) if f > 1]
            
            if keep and len(keep) < len(self.NGS_fluxes_asterism):
                self.currentAsterismIndices = [self.currentAsterismIndices[j] for j in keep]
                self.setAsterismData()
                
            self.LO_tiptilt_cov = _host(self.mLO.computeTotalResidualMatrixI(
                self.currentAsterismIndices, self.cartSciencePointingCoords,
                np.asarray(self.cartNGSCoords_asterism), self.NGS_fluxes_asterism))
            
        self.Ctot = self.LO_tiptilt_cov
        self.LO_res = np.sqrt(np.maximum(
            np.trace(self.Ctot, axis1=-2, axis2=-1), 0))
        
        if not self.addFocusError:
            return
        
        if astIndex is None:
            self.focus_cov = _host(self.mLO.computeFocusTotalResidualMatrix(
                self.cartNGSCoords_field, self.Focus_fluxes_field,
                self.Focus_freqs_field, self.Focus_SR_field,
                self.Focus_EE_field, self.Focus_FWHM_mas_field))
            
            self.GFinPSD = True
        else:
            self.focus_cov = _host(self.mLO.computeFocusTotalResidualMatrixI(
                self.currentAsterismIndices,
                np.asarray(self.cartNGSCoords_asterism), self.Focus_fluxes_asterism))
            
            self.GFinPSD = False
            
        self.GF_res = float(np.sqrt(max(self.focus_cov[0], 0)))
        
        if self.GFinPSD:
            self.PSD = self.PSD + self.GF_res**2 * self._focus_filter()
            self.HO_res = self.PSD.sum(dim=(-2, -1)).clamp_min(0).sqrt().cpu().numpy()[:, 0]


    def _combined_jitter(self):
        count = self.nPointings
        zero = torch.zeros(count, device=self.model.device, dtype=self.model.wvl.dtype)
        static = (
            fwhm_to_jitter(
                self.jitter_FWHM, device=self.model.device,
                dtype=self.model.wvl.dtype, count=count
            )
            if self.jitter_FWHM is not None else None
        )
        
        lo = None
        if self.LOisOn and self.LO_tiptilt_cov is not None:
            covariance = torch.as_tensor(self.LO_tiptilt_cov, device=self.model.device, dtype=self.model.wvl.dtype)
            lo = tiptilt_covariance_to_jitter(covariance, 2 * self.tel_radius)
            self.cov_ellipses = np.column_stack([lo[2].detach().cpu().numpy(), lo[0].detach().cpu().numpy(), lo[1].detach().cpu().numpy()])
            
        if self.LOisOn and not self.doConvolve:
            return zero, zero, zero
        
        if lo is None:
            return static if static is not None else (zero, zero, zero)
        
        if static is None:
            return lo
        
        return combine_zero_centered_jitters(*lo, *static)


    def finalPSF(self, astIndex=None):
        jitter = self._combined_jitter()
        
        if self.LOisOn and not self.doConvolveAsterism:
            if self._ho_fwhm is None:
                zeros = torch.zeros(self.nPointings, device=self.model.device,
                                    dtype=self.model.wvl.dtype)
                self.model.Jx = self.model.Jy = self.model.Jxy = zeros
                with torch.no_grad():
                    baseline = self.model(PSD=self._PSD_HO).detach().cpu().numpy()
                self._ho_fwhm = [[
                    float(np.sqrt(np.prod(psf_fwhm(baseline[j, i], self.psInMas))))
                    for j in range(self.nPointings)] for i in range(self.nWvl
                )]
                
                if self.nWvl == 1:
                    self._ho_fwhm = self._ho_fwhm[0]
                    
            self.pointings_FWHM_mas = self._ho_fwhm
            self.results = []
            return
        
        self.model.Jx, self.model.Jy, self.model.Jxy = jitter
        
        with torch.no_grad():
            self.PSF_tensor = self.model(PSD=self.PSD).detach()
            
        array = self.PSF_tensor.cpu().numpy()
        self.cubeResultsArray = (np.moveaxis(array, 1, 0) if self.nWvl > 1 else array[:, 0])
        self.cubeResults = self.cubeResultsArray
        self.results = self.cubeResultsArray
        
        self.pointings_FWHM_mas = [[
            float(np.sqrt(np.prod(psf_fwhm(array[j, i], self.psInMas))))
            for j in range(self.nPointings)] for i in range(self.nWvl)]
        
        if self.nWvl == 1:
            self.pointings_FWHM_mas = self.pointings_FWHM_mas[0]
            
        self.final_Jx_mas  = jitter[0].detach().cpu().numpy()
        self.final_Jy_mas  = jitter[1].detach().cpu().numpy()
        self.final_Jxy_deg = jitter[2].detach().cpu().numpy()


    def doOverallSimulation(self, astIndex=None):
        if self.LOisOn:
            self.configLO(astIndex)
        
        self.results = []
        
        if astIndex is None or self.firstSimCall:
            self._prepare_state()
            
        self.PSD = self._PSD_HO.clone()
        self.HO_res = self.PSD.sum(dim=(-2, -1)).clamp_min(0).sqrt().cpu().numpy()[:, 0]
        self._compute_LO(astIndex)
        self.finalPSF(astIndex)
        self.firstSimCall = False
        
        if astIndex is None and len(self.results):
            with torch.no_grad():
                self.psfOL = self.model.OLPSF(
                    include_static=True, include_jitter=False).detach().cpu().numpy()
                self.psfDL = self.model.DLPSF().detach().cpu().numpy()
            self.computePSF1D()
            
        if self.doPlot and len(self.results):
            import matplotlib.pyplot as plt
            plt.imshow(self.cubeResultsArray[0] if self.nWvl == 1
                       else self.cubeResultsArray[0, 0])
            plt.show()


    def computeMetrics(self):
        self.penalty = []
        for i, ho in enumerate(self.HO_res):
            variance = ho**2
            if self.LOisOn:
                variance += self.LO_res[i]**2
                
            if self.addFocusError and not self.GFinPSD:
                
                variance += self.GF_res**2
            self.penalty.append(float(np.sqrt(variance)))
            
        if not len(self.results):
            self.sr = [[float(np.exp(
                -4 * np.pi**2 * penalty**2 / (wavelength * 1e9)**2))
                for penalty in self.penalty] for wavelength in self.wvl]
            fwhm_LO = np.zeros(self.nPointings)
            if self.LOisOn:
                covariance = torch.as_tensor(self.LO_tiptilt_cov, dtype=torch.float64)
                jx, jy, _ = tiptilt_covariance_to_jitter(covariance, 2 * self.tel_radius)
                fwhm_LO = 2.354820045 * torch.sqrt(jx * jy).numpy()
                
            self.fwhm = [
                np.sqrt(np.asarray(widths)**2 + fwhm_LO**2).tolist()
                for widths in (self.pointings_FWHM_mas if self.nWvl > 1 else [self.pointings_FWHM_mas])
            ]
            
            self.ee = [[0.0] * self.nPointings for _ in self.wvl]
            if self.nWvl == 1:
                self.sr, self.fwhm, self.ee = self.sr[0], self.fwhm[0], self.ee[0]
            return
        
        sr, widths, energy = [], [], []
        dl = np.asarray(self.psfDL)
        
        for wave_index in range(self.nWvl):
            cube = (self.cubeResultsArray[wave_index]
                    if self.nWvl > 1 else self.cubeResultsArray)
            reference = dl[:, wave_index]
            sr_wave, width_wave, ee_wave = [], [], []
            for i, image in enumerate(cube):
                sr_wave.append(float(image.max() / reference[i].max()
                                     * reference[i].sum() / image.sum()))
                width_wave.append(float(np.sqrt(np.prod(psf_fwhm(image, self.psInMas)))))
                
                if self.ensquaredEnergy:
                    curve = psf_ensquared_energy(image)
                    radii = (np.arange(len(curve)) + 0.5) * self.psInMas
                else:
                    curve, radii = psf_encircled_energy(image, self.psInMas)
                    
                ee_wave.append(float(np.interp(self.eeRadiusInMas, radii, curve)))
                
            sr.append(sr_wave)
            widths.append(width_wave)
            energy.append(ee_wave)
        
        self.sr = sr if self.nWvl > 1 else sr[0]
        self.fwhm = widths if self.nWvl > 1 else widths[0]
        self.ee = energy if self.nWvl > 1 else energy[0]


    def computePSF1D(self):
        cube = (self.cubeResultsArray if self.nWvl > 1 else self.cubeResultsArray[None])
        y, x = np.indices(cube.shape[-2:])
        center = (cube.shape[-1] - 1) / 2
        bins = np.floor(np.hypot(x - center, y - center)).astype(int)
        count = np.bincount(bins.ravel())
        
        profiles = []
        for images in cube:
            per_wave = []
            for image in images:
                radial = np.bincount(bins.ravel(), weights=image.ravel()) / count
                per_wave.append(radial)
            profiles.append(per_wave)
            
        self.psf1d = np.asarray(profiles)
        self.psf1d_radius = np.arange(self.psf1d.shape[-1]) * self.psInMas
        self.psf1d_data = np.concatenate((np.broadcast_to(self.psf1d_radius, self.psf1d.shape), self.psf1d), axis=0)


    def savePSFprofileJSON(self):
        output = Path(self.outputDir) / f'{self.outputFile}1D_PSF.json'
        output.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            'execution_infos': {'TIME': datetime.now().strftime('%Y%m%d_%H%M%S')},
            'infos': self.my_data_map,
            'psf': {
                'radius': self.psf1d_radius.tolist(),
                'psf': self.psf1d.tolist()
            },
        }
        output.write_text(json.dumps(payload, default=str), encoding='utf-8')

    def saveResults(self):
        if not hasattr(self, 'psfOL'):
            raise RuntimeError('Run doOverallSimulation() without astIndex before saving')
        
        hdus = [
            fits.PrimaryHDU(),
            fits.ImageHDU(self.cubeResultsArray, name='PSF_CUBE'),
            fits.ImageHDU(self.psfOL, name='OPEN_LOOP'),
            fits.ImageHDU(self.psfDL, name='DIFFRACTION_LIMITED')
        ]
        
        if self.savePSDs:
            hdus.append(fits.ImageHDU(
                self.PSD.detach().cpu().numpy(), name='HIGH_ORDER_PSD'))
            
        hdus.append(fits.ImageHDU(self.psf1d_data, name='PSF_PROFILES'))
        header = hdus[1].header
        header['PIX_MAS'] = self.psInMas
        header['NSRC']    = self.nPointings
        header['NWVL']    = self.nWvl
        
        for i, wavelength in enumerate(self.wvl):
            header[f'WL_NM{i:03d}'] = wavelength * 1e9
            
        for i, residual in enumerate(self.HO_res):
            header[f'RESH{i:04d}'] = float(residual)
            
        if self.LOisOn:
            for i, residual in enumerate(self.LO_res):
                header[f'RESL{i:04d}'] = float(residual)
                
        if self.addFocusError:
            header['RESF0000'] = float(self.GF_res)
            
        if self.addSrAndFwhm:
            self.computeMetrics()
            for wave_index in range(self.nWvl):
                widths = self.fwhm[wave_index] if self.nWvl > 1 else self.fwhm
                strehls = self.sr[wave_index] if self.nWvl > 1 else self.sr
                for i, (width, strehl) in enumerate(zip(widths, strehls)):
                    header[f'FW{wave_index:02d}{i:03d}'] = float(width)
                    header[f'SR{wave_index:02d}{i:03d}'] = float(strehl)
                    
        output = Path(self.outputDir) / f'{self.outputFile}.fits'
        output.parent.mkdir(parents=True, exist_ok=True)
        fits.HDUList(hdus).writeto(output, overwrite=True)
