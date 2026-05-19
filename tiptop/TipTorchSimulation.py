# import p3.aoSystem
# import p3.aoSystem.fourierModel
# import p3.aoSystem.FourierUtils
# from p3.aoSystem.fourierModel import *
# from p3.aoSystem.FourierUtils import *

import os
import copy
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"


from tiptorch.PSF_models.TipTorch import TipTorch
from tiptorch.managers.config_manager import ConfigManager
from tiptorch._config import default_device, default_torch_type
import torch

try:
    from utils import cov_to_jitter_mas
except Exception:
    cov_to_jitter_mas = None

try:
    # In the user project this helper may live in RotateGaussian.py.
    from RotateGaussian import combine_zero_centered_jitters
except Exception:
    try:
        from Gauss_rotate import combine_zero_centered_jitters
    except Exception:
        combine_zero_centered_jitters = None

CALIBRATIONS_PATH = os.path.normpath(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 
        '../P3/p3/aoSystem/data/'
    )
)

# import json

from mastsel import *

# from .tiptopUtils import *
# from ._version import __version__
# import matplotlib as mpl


rad2mas = 3600 * 180 * 1000 / np.pi

class baseSimulation(object):

    def __init__(
        self,
        path: str,
        parametersFile: str,
        outputDir: str,
        outputFile: str,
        doConvolve: bool = True,
        doPlot: bool = False,
        addSrAndFwhm: bool = True,
        verbose: bool = False,
        getHoErrorBreakDown: bool = False,
        savePSDs: bool = False,
        ensquaredEnergy: bool = False,
        eeRadiusInMas: float = 50.0
    ):

        self.path = path
        self.verbose = verbose
        self.firstSimCall = True    
        self.doConvolveAsterism = True
        self.pointings_FWHM_mas = None
        self.parametersFile = parametersFile
        self.outputDir = outputDir
        self.outputFile = outputFile
        self.doPlot = doPlot
        self.doConvolve = doConvolve
        self.addSrAndFwhm = addSrAndFwhm
        self.getHoErrorBreakDown = getHoErrorBreakDown
        self.savePSDs = savePSDs
        self.ensquaredEnergy = ensquaredEnergy
        self.eeRadiusInMas = eeRadiusInMas
        
        if verbose:
            np.set_printoptions(precision=3)
 
        config_manager = ConfigManager()
        config_dict = config_manager.Load( os.path.normpath(os.path.join(self.path, self.parametersFile)) )
        config_dict = config_manager.Convert(config_dict, framework='pytorch', device=default_device, dtype=default_torch_type)
        self.config_torch = config_dict
        
        # Initialize TipTorch model
        self.model = TipTorch(
            AO_config = config_dict,
            norm_regime = 'sum',
            device = default_device,
            oversampling = 1,
            retain_PSDs = True
        )
        
        self.overSamp = self.model.oversampling # Oversampling factor
        self.psInMas = float(self.model.psInMas.detach().cpu().item() if torch.is_tensor(self.model.psInMas) else self.model.psInMas)
        self.nPixPSF = int(self.model.N_pix.detach().cpu().item() if torch.is_tensor(self.model.N_pix) else self.model.N_pix)
        self.nPointings = int(self.model.N_src)
        self.wvl = self.model.wvl.flatten().cpu().numpy().tolist()
        self.wvlMax = max(self.wvl)
        self.nWvl = len(self.wvl)
        
        # Convert back from torch
        self.config = config_manager.Convert(config_dict, framework='list')
        
        self.tel_radius = self.config['telescope']['TelescopeDiameter'] / 2.0  # [m]
            
        self.zenithSrc  = self.config['sources_science']['Zenith']
        self.azimuthSrc = self.config['sources_science']['Azimuth']
        self.pointings  = polarToCartesian(np.array( [self.zenithSrc, self.azimuthSrc]))
        
        self.xxSciencePointigs = self.pointings[0,:]
        self.yySciencePointigs = self.pointings[1,:]
        
        # self.psInMas = self.config['sensor_science']['PixelScale']
        # self.SupSamp = self.config['sensor_science'].get('Super_Sampling', None)
        
        # it checks if LO parameters are set and then it acts accordingly
        if 'sensor_LO' in self.config.keys():
            self.LOisOn = True
            if self.verbose:
                print('LO part is present')
        else:
            self.LOisOn = False
            self.nNaturalGS_field = 0
            
            if self.verbose:
                print('LO part is not present')
    
        self.jitter_FWHM   = self.config['telescope'].get('jitter_FWHM', None)
        self.addFocusError = self.config['telescope'].get('glFocusOnNGS', False)
        
        self.GFinPSD = False # defines whether the optional focus error is added as a fixed term in the PSD (True) or as a post-PSD quadratic term (False)

        # LO covariance products are computed later, per simulation call.
        # LO_tiptilt_cov: tip/tilt-like image-motion covariance used for PSF convolution.
        # focus_cov: optional global-focus covariance reduced to GF_res.
        self.LO_tiptilt_cov = None
        self.focus_cov = None
        
        if not 'sensor_Focus' in self.config.keys() and self.addFocusError and max(self.config['sensor_LO']['NumberLenslets']) == 1:
            raise ValueError("[telescope] glFocusOnNGS (that is focus correction with NGS) is available only if NGS/Focus WFSs have more than one sub-aperture")


    def configLO(self):
        def first_if_list(value):
            return value[0] if isinstance(value, list) else value

        def as_list(value, n):
            return value if isinstance(value, list) else [value] * n

        config = self.config
        sources_LO = config["sources_LO"]
        sensor_LO  = config["sensor_LO"]
        rtc = config.get("RTC", {})

        self.cartSciencePointingCoords = np.dstack( (self.xxSciencePointigs, self.yySciencePointigs) ).reshape(-1, 2)

        self.LO_zen_field = sources_LO["Zenith"]
        self.LO_az_field  = sources_LO["Azimuth"]

        N_LO = len(self.LO_zen_field)

        self.LO_wvl          = first_if_list(sources_LO["Wavelength"])
        self.LO_fluxes_field = sensor_LO["NumberPhotons"]
        self.LO_psInMas      = as_list(sensor_LO["PixelScale"], N_LO)
        self.LO_freqs_field  = as_list(rtc["SensorFrameRate_LO"], N_LO)

        self.addLOAlias = sensor_LO.get("addAliasError", False)

        sensor_focus = config.get("sensor_Focus")

        if sensor_focus is not None:
            sources_focus = config.get("sources_Focus", sources_LO)

            self.Focus_fluxes4s_field = sensor_focus["NumberPhotons"]
            self.Focus_psInMas = as_list(sensor_focus["PixelScale"], N_LO)
            self.Focus_wvl = first_if_list(sources_focus["Wavelength"])
        else:
            self.Focus_fluxes4s_field = self.LO_fluxes_field
            self.Focus_psInMas = self.LO_psInMas
            self.Focus_wvl = self.LO_wvl

        self.Focus_freqs_field = as_list( rtc.get("SensorFrameRate_Focus", self.LO_freqs_field), N_LO, )

        self.NGS_fluxes_field   = [flux*freq for flux, freq in zip(self.LO_fluxes_field, self.LO_freqs_field)]
        self.Focus_fluxes_field = [flux*freq for flux, freq in zip(self.Focus_fluxes4s_field, self.Focus_freqs_field)]

        polarNGSCoords = np.column_stack((self.LO_zen_field, self.LO_az_field))

        self.nNaturalGS_field = N_LO
        self.cartNGSCoords_field = np.asarray([ polarToCartesian(coords) for coords in polarNGSCoords ])
        self.currentAsterismIndices = list(range(N_LO))

        indices = self.currentAsterismIndices

        self.LO_zen_asterism        = [self.LO_zen_field[i] for i in indices]
        self.LO_az_asterism         = [self.LO_az_field[i] for i in indices]
        self.LO_fluxes_asterism     = [self.LO_fluxes_field[i] for i in indices]
        self.LO_freqs_asterism      = [self.LO_freqs_field[i] for i in indices]
        self.NGS_fluxes_asterism    = [self.NGS_fluxes_field[i] for i in indices]
        self.Focus_fluxes_asterism  = [self.Focus_fluxes_field[i] for i in indices]
        self.cartNGSCoords_asterism = [self.cartNGSCoords_field[i] for i in indices]
    
    
    '''
    def finalPSF(self,astIndex):
        if astIndex is None or self.firstSimCall:
            # ----------------------------------------------------------------------------
            ## HO PSF
            PSD_HO = arrayP3toMastsel(self.PSD[0:self.nPointings])
            mask = arrayP3toMastsel(self.fao.ao.tel.pupil)
            padPSD = self.nWvl > 1

            if self.verbose:
                print('******** HO PSF')
                
            psfLongExpPointingsArr = psdSetToPsfSet(PSD_HO, mask,
                                                    self.wvl, self.N, self.sx, self.grid_diameter,
                                                    self.freq_range, self.dk, self.nPixPSF,
                                                    self.wvlMax, self.overSamp,
                                                    opdMap=self.opdMap, padPSD=padPSD)
            # -----------------------------------------------------------------
            ## Merit functions
            self.pointings_FWHM_mas   = []

            for i in range(self.nWvl):
                if self.nWvl>1:
                    psfList = psfLongExpPointingsArr[i]
                    wvl = self.wvl[i]
                else:
                    psfList = psfLongExpPointingsArr
                    wvl = self.wvl[0]
                fwhmList = []
                idx = 0
                for img in psfList:
                    # Get SFWHM in mas the star positions at the sensing wavelength
                    fwhmX,fwhmY = getFWHM(img.sampling, self.psInMas, method='contour', nargout=2)
                    fwhm = np.sqrt(fwhmX*fwhmY)
                    fwhmList.append(fwhm) #average over major and minor axes
                    if self.verbose:
                        s1 = cpuArray(PSD_HO[idx]).sum()
                        sr = np.exp(-s1*(2*np.pi*1e-9/wvl)**2) # Strehl-ratio at the sensing wavelength
                        print('SR(@',int(wvl*1e9),'nm)        :', "%.5f" % sr)
                        print('FWHM(@',int(wvl*1e9),'nm) [mas]:', "%.3f" % fwhm)
                    idx += 1
                if self.nWvl > 1:
                    self.pointings_FWHM_mas.append(fwhmList)
                else:
                    self.pointings_FWHM_mas = fwhmList

            self.psfLongExpPointingsArr = psfLongExpPointingsArr

            # ----------------------------------------------------------------------------
            ## computation of the HO error (this is fixed for the simulation)
            self.HO_res = np.sqrt(np.sum(self.PSD[0:self.nPointings],axis=(1,2)))

        def jitter_ellipse():
            if self.jitter_FWHM is None:
                return None
            
            if isinstance(self.jitter_FWHM, list):
                return np.array([
                    self.jitter_FWHM[2],
                    sigma_from_FWHM(self.jitter_FWHM[0]),
                    sigma_from_FWHM(self.jitter_FWHM[1])
                ], dtype=self.fao.dtype)
            
            return np.array([
                0,
                sigma_from_FWHM(self.jitter_FWHM),
                sigma_from_FWHM(self.jitter_FWHM)
            ], dtype=self.fao.dtype)


        def ellipse_to_spectrum(ellp):
            if ellp is None:
                return None
            return residualToSpectrum(
                ellp, self.wvlRef, self.nPixPSF,
                1/(self.nPixPSF * self.psInMas)
            )


        def store_results(convolve_one):
            for i in range(self.nWvl):
                psfList = self.psfLongExpPointingsArr[i] if self.nWvl > 1 else self.psfLongExpPointingsArr
                resultList = [convolve_one(psf, idx) for idx, psf in enumerate(psfList)]
                if self.nWvl > 1:
                    self.results.append(resultList)
                else:
                    self.results = resultList

        # ------------------------------------------------------------------------
        ## final PSFs computation after optional convolution with LO/jitter kernels
        jitterSpec = ellipse_to_spectrum(jitter_ellipse())

        if self.LOisOn:
            if not self.doConvolve:
                store_results(lambda psf, idx: psf)
                return

            self.cov_ellipses = self.mLO.ellipsesFromCovMats(self.LO_tiptilt_cov)

            if self.verbose:
                for n in range(self.cov_ellipses.shape[0]):
                    print('cov_ellipses #', n, ': ', self.cov_ellipses[n,:], ' (unit: rad, mas, mas)')

            # Preserve old behavior: cache ellipses only when asterism convolution is disabled.
            if not self.doConvolveAsterism:
                return

            if self.verbose:
                print('******** FINAL CONVOLUTION')

            loSpecs = [
                ellipse_to_spectrum(ellp.astype(self.fao.dtype))
                for ellp in self.cov_ellipses
            ]

            def convolve_one(psf, idx):
                psf = convolve(psf, loSpecs[idx])
                return convolve(psf, jitterSpec) if jitterSpec is not None else psf

            store_results(convolve_one)
        else:
            store_results(
                lambda psf, idx: convolve(psf, jitterSpec) if jitterSpec is not None else psf
            )
    '''
    
    
    def _sensor_PSF_sampling(self, sensor_wvl, sensor_pixel_scales):
        """Return PSF pixel scale, output size, and reshape mode for an LO-like sensor.

        The NGS/focus PSFs are generated at the science-camera sampling scaled to
        the sensor wavelength.  If that sampling is coarser than the sensor pixel
        scale, keep the native high-resolution PSF instead of reshaping too early.
        """
        psf_pixel_scale_mas = self.model.psInMas * sensor_wvl / self.wvlMax
        min_sensor_pixel_scale = np.min(sensor_pixel_scales)

        if psf_pixel_scale_mas / min_sensor_pixel_scale > 1 and self.overSamp > 1:
            return psf_pixel_scale_mas / self.overSamp, int(self.overSamp * self.nPixPSF), True

        return psf_pixel_scale_mas, self.nPixPSF, False


    def _build_sensor_PSD_and_mask(self, n_lenslets):
        """Build NGS-direction PSDs and the matching LO/focus sub-aperture mask."""
        PSD   = arrayP3toMastsel(self.PSD[-self.nNaturalGS_field:])
        pupil = arrayP3toMastsel(self.fao.ao.tel.pupil)
        mask  = maskSA(n_lenslets, self.nNaturalGS_field, pupil)
        self._apply_subaperture_piston_filter(PSD, n_lenslets)
        
        return PSD, mask


    def _apply_subaperture_piston_filter(self, PSD, n_lenslets):
        """Filter each NGS PSD by the piston filter of its effective sub-aperture.

        A one-lenslet sensor is treated as a pure full-aperture measurement and is
        left unchanged.  If the config provides one lenslet count for all stars,
        that value is reused for every NGS direction.
        """
        spatial_freq = np.sqrt(self.fao.freq.k2_)

        for i in range(self.nNaturalGS_field):
            n_lenslets_i = n_lenslets[i] if len(n_lenslets) == self.nNaturalGS_field else n_lenslets[0]

            if n_lenslets_i == 1:
                continue

            piston_filter = FourierUtils.pistonFilter(2 * self.tel_radius / n_lenslets_i, spatial_freq)
            PSD[i] *= piston_filter


    def _PSFs_from_sensor_PSDs(self, PSD, mask, sensor_wvl, nPixPSF, skip_reshape):
        """Convert a set of NGS-direction PSDs into long-exposure sensor PSFs."""
        return psdSetToPsfSet(
            PSD,
            mask,
            sensor_wvl,
            self.N,
            self.sx,
            self.grid_diameter,
            self.freq_range,
            self.dk,
            nPixPSF,
            self.wvlMax,
            self.overSamp,
            opdMap=self.opdMap,
            skip_reshape=skip_reshape
        )


    def _encircled_energy_at_sensor_radius(self, img, fwhm_mas, psf_pixel_scale_mas, sensor_pixel_scale_mas, nPixPSF):
        """Measure EE at max(FWHM, sensor pixel scale), matching the old logic."""
        if 2 * fwhm_mas >= nPixPSF * psf_pixel_scale_mas:
            return 1

        ee, rr = getEncircledEnergy(
            img.sampling,
            pixelscale=psf_pixel_scale_mas,
            center=centeredPixelCoords(nPixPSF),
            nargout=2
        )
        ee *= 1 / np.max(ee)
        ee_at_radius = interp1d(rr, ee, kind='cubic', bounds_error=False)
        return ee_at_radius(max([fwhm_mas, sensor_pixel_scale_mas]))


    def _measure_sensor_PSFs(self, PSFs, PSD, sensor_wvl, PSF_pixel_scale_mas, sensor_pixel_scales, nPixPSF, label):
        """Return SR, FWHM, and EE lists for NGS/focus sensor PSFs."""
        sr_list, fwhm_list, ee_list = [], [], []

        for idx, img in enumerate(PSFs):
            psd_variance_nm2 = cpuArray(PSD[idx]).sum()
            SR = np.exp(-psd_variance_nm2 * (2 * np.pi * 1e-9 / sensor_wvl)**2)

            fwhmX, fwhmY = getFWHM(img.sampling, PSF_pixel_scale_mas, method='contour', nargout=2)
            FWHM = np.sqrt(fwhmX * fwhmY)

            sensor_pixel_scale_i = sensor_pixel_scales[idx] if isinstance(sensor_pixel_scales, list) else sensor_pixel_scales
            EE = self._encircled_energy_at_sensor_radius(
                img,
                FWHM,
                PSF_pixel_scale_mas,
                sensor_pixel_scale_i,
                nPixPSF
            )

            sr_list.append(SR)
            fwhm_list.append(FWHM)
            ee_list.append(EE)

            if self.verbose:
                print('SR(@', int(sensor_wvl * 1e9), 'nm)        :', "%.5f" % SR)
                print('FWHM(@', int(sensor_wvl * 1e9), 'nm) [mas]:', "%.3f" % FWHM)
                print(label, ':', "%.5f" % EE)

        return sr_list, fwhm_list, ee_list


    def _compute_NGS_DL_FWHM(self, maskLO, n_lenslet_masks, psf_pixel_scale_mas):
        """Compute diffraction-limited NGS FWHM used by the optional LO aliasing term."""
        if not self.addLOAlias:
            self.NGS_DL_FWHM_mas = None
            return

        if self.verbose:
            print('Adding aliasing error on LO!')

        self.NGS_DL_FWHM_mas = [] if isinstance(maskLO, list) else None

        for i in range(n_lenslet_masks):
            mask_i = maskLO[i] if isinstance(maskLO, list) else maskLO

            PSD_DL    = Field(self.LO_wvl, self.N, self.freq_range, 'rad')
            maskField = Field(self.LO_wvl, self.N, self.grid_diameter)
            
            maskField.sampling = congrid(mask_i, [self.sx, self.sx])
            maskField.sampling = padOrCropCentered(maskField.sampling, self.N, xp=maskField.xp)

            psfNgsDL = longExposurePsf(maskField, PSD_DL)
            fwhmX, fwhmY = getFWHM(psfNgsDL.sampling, psf_pixel_scale_mas, method='contour', nargout=2)
            fwhm = np.sqrt(fwhmX * fwhmY)

            if self.NGS_DL_FWHM_mas is None:
                self.NGS_DL_FWHM_mas = fwhm
            else:
                self.NGS_DL_FWHM_mas.append(fwhm)


    def _use_NGS_metrics_for_focus(self):
        """Reuse LO NGS PSF metrics when no separate focus sensor is configured."""
        if self.verbose:
            print('Focus sensor is not set: using LO PSFs.')

        self.Focus_SR_field = self.NGS_SR_field
        self.Focus_FWHM_mas_field = self.NGS_FWHM_mas_field
        self.Focus_EE_field = self.NGS_EE_field


    def _compute_focus_sensor_PSF_metrics(self):
        """Compute PSF quality metrics for a dedicated focus sensor, if present."""
        if not self.addFocusError:
            return

        if 'sensor_Focus' not in self.config:
            self._use_NGS_metrics_for_focus()
            return

        if self.verbose:
            print('Focus sensor is set: computing new PSFs.')
            print('******** Focus Sensor PSF - NGS directions (1 sub-aperture)')

        PSF_pixel_scale_mas, nPixPSF, skip_reshape = self._sensor_PSF_sampling(self.Focus_wvl, self.Focus_psInMas)

        n_lenslets = self.config['sensor_Focus']['NumberLenslets']
        PSD_focus, mask_focus = self._build_sensor_PSD_and_mask(n_lenslets)
        
        PSF_focus = self._PSFs_from_sensor_PSDs(
            PSD_focus,
            mask_focus,
            self.Focus_wvl,
            nPixPSF,
            skip_reshape
        )

        self.Focus_SR_field, self.Focus_FWHM_mas_field, self.Focus_EE_field = self._measure_sensor_PSFs(
            PSF_focus,
            PSD_focus,
            self.Focus_wvl,
            PSF_pixel_scale_mas,
            self.Focus_psInMas,
            nPixPSF,
            label='EE (focus sensor)    '
        )


    def ngsPSF(self):
        """Compute PSF-based quantities needed by the low-order MavisLO model.

        This does not produce the final science PSF.  It prepares guide-star
        quality inputs for the LO covariance calculation:
        1. Generate long-exposure PSFs in the NGS directions.
        2. Measure NGS SR/FWHM/EE from those PSFs.
        3. Optionally compute diffraction-limited NGS FWHM for LO aliasing.
        4. If focus correction uses a dedicated sensor, repeat the same metric
           calculation for that sensor; otherwise reuse the NGS metrics.
        """
        if self.verbose:
            print('******** LO PSF - NGS directions (1 sub-aperture)')

        PSF_pixel_scale_mas, nPixPSF, skip_reshape = self._sensor_PSF_sampling(self.LO_wvl, self.LO_psInMas)

        n_lenslets = self.config['sensor_LO']['NumberLenslets']
        n_lenslet_masks = len(n_lenslets) if len(n_lenslets) == self.nNaturalGS_field else 1

        PSD_NGS, mask_LO = self._build_sensor_PSD_and_mask(n_lenslets)
        PSF_NGS = self._PSFs_from_sensor_PSDs(PSD_NGS, mask_LO, self.LO_wvl, nPixPSF, skip_reshape)

        self.NGS_SR_field, self.NGS_FWHM_mas_field, self.NGS_EE_field = self._measure_sensor_PSFs(
            PSF_NGS,
            PSD_NGS,
            self.LO_wvl,
            PSF_pixel_scale_mas,
            self.LO_psInMas,
            nPixPSF,
            label = 'EE                  '
        )

        self._compute_NGS_DL_FWHM(mask_LO, n_lenslet_masks, PSF_pixel_scale_mas)
        self._compute_focus_sensor_PSF_metrics()


    '''
    def computeMetrics(self):
        self.penalty, self.sr, self.fwhm, self.ee = [], [], [], []
        
        if len(self.results) == 0:
            if self.LOisOn:
                if self.addFocusError and not self.GFinPSD:
                    self.penalty.append( np.sqrt( np.mean(cpuArray(self.LO_res)**2 + cpuArray(self.HO_res)**2) + self.GF_res**2 ) )
                else:
                    self.penalty.append( np.sqrt( np.mean(cpuArray(self.LO_res)**2 + cpuArray(self.HO_res)**2) ) )
            else:
                self.penalty.append( np.sqrt( np.mean(cpuArray(self.HO_res)**2) ) )
                
            self.sr.append( np.exp( -4*np.pi**2 * ( self.penalty[-1]**2 )/(self.wvlRef*1e9)**2) )
            scale = (np.pi/(180*3600*1000) * 2 * self.tel_radius / (4*1e-9))
            FWHMS_LO = 2.355 * self.LO_res/scale / np.sqrt(2)
            self.fwhm.append(np.sqrt(FWHMS_LO**2 + np.asarray(self.pointings_FWHM_mas)**2))
            self.ee.append(0)
            
        else:
            for idx, HO_res in enumerate(cpuArray(self.HO_res)):
                if self.LOisOn:
                    if self.addFocusError and not self.GFinPSD:
                        self.penalty.append( np.sqrt( cpuArray(self.LO_res)[idx]**2 + HO_res**2 + self.GF_res**2) )
                    else:
                        self.penalty.append( np.sqrt( cpuArray(self.LO_res)[idx]**2 + HO_res**2 ) )
                else:
                    self.penalty.append( HO_res )
            if self.verbose:
                print('EE is computed for a radius of ', self.eeRadiusInMas,' mas')

            for i in range(self.nWvl):
                if self.nWvl > 1:
                    results = self.results[i]
                else:
                    results = self.results
                    
                samp = self.wvl[i] * rad2mas / (self.psInMas*2*self.tel_radius)
                
                sr, fwhm, ee = [], [], []
                
                for img in results:
                    sr.append(getStrehl(img.sampling, self.fao.ao.tel.pupil, samp, method='max', psfInOnePix=True))
                    fwhm.append(getFWHM(img.sampling, self.psInMas, method='contour', nargout=1))
                    
                    if self.ensquaredEnergy:
                        ee_ = cpuArray(getEnsquaredEnergy(img.sampling))
                        rr_ = np.arange(1, ee_.shape[0]*2, 2) * self.psInMas * 0.5
                    else:
                        ee_,rr_ = getEncircledEnergy(img.sampling, pixelscale=self.psInMas,
                                                     center=centeredPixelCoords(self.nPixPSF), nargout=2)
                    ee_at_radius_fn = interp1d(rr_, ee_, kind='cubic', bounds_error=False)
                    ee.append( cpuArray(ee_at_radius_fn(self.eeRadiusInMas)).item() )
                    
                if self.nWvl > 1:
                    self.sr.append(sr)
                    self.fwhm.append(fwhm)
                    self.ee.append(ee)
                else:
                    self.sr = sr
                    self.fwhm = fwhm
                    self.ee = ee
    '''

    # def _configure_P3_LO_maps(self):
    #     """Push LO source/sensor values into the P3/fourierModel structures."""
    #     if 'sensor_LO' in self.config:
    #         photons = self.config['sensor_LO']['NumberPhotons']
    #         self.fao.my_data_map['sensor_LO']['NumberPhotons'] = photons
    #         self.fao.ao.my_data_map['sensor_LO']['NumberPhotons'] = photons

    #     if 'sources_LO' in self.config:
    #         self.fao.my_data_map['sources_LO'] = self.config['sources_LO']
    #         self.fao.ao.my_data_map['sources_LO'] = self.config['sources_LO']
    #         self.fao.ao.configLOsensor()
    #         self.fao.ao.configLO()
    #         self.fao.ao.configLO_SC()
    

    # def _init_HO_fourier_model(self):
    #     """Initialize P3/fourierModel and cache the geometry needed downstream."""
    #     if self.verbose:
    #         print('******** HO PSD science and NGSs directions')

    #     self.fao = fourierModel(
    #         self.fullPathFilename,
    #         calcPSF=False,
    #         verbose=self.verbose,
    #         display=False,
    #         getPSDatNGSpositions=self.LOisOn,
    #         computeFocalAnisoCov=False,
    #         TiltFilter=self.LOisOn,
    #         getErrorBreakDown=self.getHoErrorBreakDown,
    #         doComputations=False,
    #         psdExpansion=True,
    #         reduce_memory=True
    #     )

    #     self._configure_P3_LO_maps()

    #     if self.verbose:
    #         print('Setting MASTSEL PSF precision to:', self.fao.dtype)

    #     mastselPsfPrecision(dtype=self.fao.dtype)
    #     self.fao.initComputations()

    #     # High-order PSD calculations at science directions and, optionally, NGS directions.
    #     self.PSD = self.fao.PSD.transpose()  # [nm^2]
    #     self.N = self.PSD[0].shape[0]
    #     self.nPointings = self.pointings.shape[1]
    #     self.nPixPSF = int(self.fao.ao.cam.fovInPix)
    #     self.overSamp = int(self.fao.freq.kRef_)
    #     self.PSDstep = self.fao.freq.PSDstep
    #     self.freq_range = self.N * self.PSDstep
    #     self.grid_diameter = 1 / self.PSDstep
    #     self.sx = int(2 * np.round(self.tel_radius * self.freq_range))

    #     # Same as in p3.aoSystem.powerSpectrumDensity, except multiplied by 1e9 instead of 2.
    #     self.dk = 1e9 * self.fao.freq.kcMax_ / self.fao.freq.resAO

    #     # P3 reference wavelength, required to scale the OL PSD from rad to m.
    #     self.wvlRef = self.fao.freq.wvlRef

    #     self.mask = Field(self.wvlRef, self.N, self.grid_diameter)
    #     self.mask.sampling = congrid(arrayP3toMastsel(self.fao.ao.tel.pupil), [self.sx, self.sx])
    #     self.mask.sampling = padOrCropCentered(self.mask.sampling, self.N, xp=self.mask.xp)

    #     if abs(float(self.psInMas) - float(cpuArray(self.fao.freq.psInMas[0]))) > 1e-6:
    #         raise ValueError(
    #             "sensor_science.PixelScale, '{}', is different from self.fao.freq.psInMas,'{}'"
    #             .format(self.psInMas, cpuArray(self.fao.freq.psInMas[0]))
    #         )

    #     if self.fao.ao.tel.opdMap_on is not None:
    #         self.opdMap = arrayP3toMastsel(self.fao.ao.tel.opdMap_on)
    #     else:
    #         self.opdMap = None

    #     if self.verbose:
    #         print('PSD step:', self.PSDstep)
    #         print('PSD freq range:', self.freq_range)
    #         print('PSD shape:', self.PSD.shape)
    #         print('PSD dtype:', self.PSD.dtype)
    #         print('oversampling:', self.overSamp)
    #         print('sensor_science.PixelScale:', self.psInMas)


    def _prepare_static_PSF_state(self, astIndex):
        """Prepare the expensive static state used by finalPSF and LO covariance calls."""
        """Return True when the static HO/NGS/MavisLO state must be recomputed."""
        if not (astIndex is None or self.firstSimCall):
            return       

        # self._init_HO_fourier_model()

        if not self.LOisOn:
            return

        if self.verbose:
            print('******** LO PART')

        # self.ngsPSF()
        self.mLO = MavisLO(self.path, self.parametersFile, verbose=self.verbose)


    def _compute_full_field_LO_terms(self):
        """Compute full-field LO terms using all available guide stars.

        Produces two separate quantities:
        - self.LO_tiptilt_cov: residual image-motion covariance in science directions.
          finalPSF converts this into anisotropic convolution ellipses/kernels.
        - self.focus_cov/self.GF_res: optional global-focus residual. In full-field
          mode it is injected into self.PSD immediately, so it is not added again
          in the residual metric budget.
        """
        self.LO_tiptilt_cov = self.mLO.computeTotalResidualMatrix(
            np.array(self.cartSciencePointingCoords),
            self.cartNGSCoords_field,
            self.NGS_fluxes_field,
            self.LO_freqs_field,
            self.NGS_SR_field,
            self.NGS_EE_field,
            self.NGS_FWHM_mas_field,
            aNGS_FWHM_DL_mas=self.NGS_DL_FWHM_mas,
            doAll=True
        )
        self.LO_res = np.sqrt(np.trace(self.LO_tiptilt_cov, axis1=1, axis2=2))

        if not self.addFocusError:
            return

        self.focus_cov = self.mLO.computeFocusTotalResidualMatrix(
            self.cartNGSCoords_field,
            self.Focus_fluxes_field,
            self.Focus_freqs_field,
            self.Focus_SR_field,
            self.Focus_EE_field,
            self.Focus_FWHM_mas_field
        )
        
        self.GF_res = np.sqrt(max(self.focus_cov[0], 0))
        self.GFinPSD = True

        # Add the global-focus residual to HO PSD
        FocusFilter = self.fao.FocusFilter()
        FocusFilter /= FocusFilter.sum()

        self.PSD += self.GF_res**2 * FocusFilter


    def _prepare_asterism_LO_cache(self):
        """Initialize MavisLO internal caches before per-asterism covariance calls."""
        if not self.firstSimCall:
            return

        # Warm up/cache the full guide-star geometry for the indexed tip/tilt calls.
        self.mLO.computeTotalResidualMatrix(
            np.array(self.cartSciencePointingCoords),
            self.cartNGSCoords_field,
            self.NGS_fluxes_field,
            self.LO_freqs_field,
            self.NGS_SR_field,
            self.NGS_EE_field,
            self.NGS_FWHM_mas_field,
            aNGS_FWHM_DL_mas=self.NGS_DL_FWHM_mas,
            doAll=False
        )

        if self.addFocusError:
            # Warm up/cache the full guide-star geometry for the indexed focus calls.
            self.mLO.computeFocusTotalResidualMatrix(
                self.cartNGSCoords_field,
                self.Focus_fluxes_field,
                self.Focus_freqs_field,
                self.Focus_SR_field,
                self.Focus_EE_field,
                self.Focus_FWHM_mas_field
            )


    def _discard_too_faint_asterism_stars(self):
        """Drop guide stars with less than one photon per frame per sub-aperture."""
        if not (np.min(self.NGS_fluxes_asterism) < 1 and np.max(self.NGS_fluxes_asterism) > 1):
            return

        valid_indices = np.where(np.array(self.NGS_fluxes_asterism) > 1)[0]
        self.NGS_fluxes_asterism    = [ elem for i, elem in enumerate(self.NGS_fluxes_asterism)    if i in valid_indices ]
        self.Focus_fluxes_asterism  = [ elem for i, elem in enumerate(self.Focus_fluxes_asterism)  if i in valid_indices ]
        self.cartNGSCoords_asterism = [ elem for i, elem in enumerate(self.cartNGSCoords_asterism) if i in valid_indices ]
        self.currentAsterismIndices = [ elem for i, elem in enumerate(self.currentAsterismIndices) if i in valid_indices ]


    def _compute_asterism_LO_terms(self):
        """Compute LO terms for the currently selected guide-star asterism.

        The tip/tilt covariance changes with the selected guide-star subset and is
        still used by finalPSF for LO convolution. The focus residual is also
        recomputed for this asterism, but it is deliberately kept outside self.PSD:
        the HO PSF is reused across asterism trials, and injecting focus into the
        cached PSD would contaminate later trials.
        """
        self._prepare_asterism_LO_cache()
        self._discard_too_faint_asterism_stars()

        self.LO_tiptilt_cov = self.mLO.computeTotalResidualMatrixI(
            self.currentAsterismIndices,
            np.array(self.cartSciencePointingCoords),
            np.array(self.cartNGSCoords_asterism),
            self.NGS_fluxes_asterism
        )
        self.LO_res = np.sqrt(np.trace(self.LO_tiptilt_cov, axis1=1, axis2=2))

        if not self.addFocusError:
            return

        self.focus_cov = self.mLO.computeFocusTotalResidualMatrixI(
            self.currentAsterismIndices,
            np.array(self.cartNGSCoords_asterism),
            self.Focus_fluxes_asterism
        )
        self.GF_res = np.sqrt(max(self.focus_cov[0], 0))
        self.GFinPSD = False


    def _compute_LO_terms(self, astIndex):
        """Compute per-call LO terms used by finalPSF and residual metrics."""
        if not self.LOisOn:
            return

        if astIndex is None:
            self._compute_full_field_LO_terms()
        else:
            self._compute_asterism_LO_terms()


    '''
    def _plot_final_PSFs(self):
        if not self.doPlot:
            return

        results = self.results[0] if self.nWvl > 1 else self.results

        if self.LOisOn and self.doConvolve:
            tiledDisplay(results)
            plotEllipses(self.cartSciencePointingCoords, self.cov_ellipses, 0.4)
        else:
            results[0].standardPlot(True)


    def _collect_cube_results(self):
        """Convert Field PSFs into plain arrays used by FITS/profile output."""
        self.cubeResults = []
        cubeResultsArray = []

        for i in range(self.nWvl):
            results = self.results[i] if self.nWvl > 1 else self.results
            cubeResults = [cpuArray(img.sampling) for img in results]

            if self.nWvl > 1:
                self.cubeResults.append(cubeResults)
                cubeResultsArray.append(np.array(cubeResults))
            else:
                self.cubeResults = cubeResults
                cubeResultsArray = cubeResults

        self.cubeResultsArray = np.array(cubeResultsArray)

    def _print_residual_summary(self):
        if not self.verbose:
            return

        print('HO_res [nm]:', self.HO_res)

        if self.LOisOn:
            print('LO_res [nm]:', self.LO_res)

        if hasattr(self, 'GF_res'):
            print('GF_res [nm]:', self.GF_res)


    def _finalize_full_field_results(self, astIndex):
        """Only full-field runs build OL/DL PSFs and radial-profile products."""
        if astIndex is not None:
            return

        self.computeOL_PSD()
        self.computeDL_PSD()
        self._collect_cube_results()
        self.computePSF1D()
        self._print_residual_summary()
    '''


    # -------------------------------------------------------------------------
    # TipTorch-native simulation path.
    # The P3-based versions above are intentionally left commented as reference.

    def _as_model_tensor(self, value, *, dtype=None):
        """Create a tensor on the same device/dtype as the TipTorch model."""
        dtype = dtype or self.model.wvl.dtype
        return torch.as_tensor(value, device=self.model.device, dtype=dtype)


    def _zero_model_jitter(self, model=None):
        """Mute TipTorch's OTF-domain jitter convolution."""
        model = self.model if model is None else model
        model.Jx  = torch.zeros(model.N_src, device=model.device, dtype=model.wvl.dtype)
        model.Jy  = torch.zeros(model.N_src, device=model.device, dtype=model.wvl.dtype)
        model.Jxy = torch.zeros(model.N_src, device=model.device, dtype=model.wvl.dtype)


    def _set_model_jitter(self, Jx, Jy, Jxy):
        """Set final jitter parameters used internally by TipTorch.JitterKernel()."""
        self.model.Jx  = Jx.to(self.model.device, dtype=self.model.wvl.dtype).flatten()
        self.model.Jy  = Jy.to(self.model.device, dtype=self.model.wvl.dtype).flatten()
        self.model.Jxy = Jxy.to(self.model.device, dtype=self.model.wvl.dtype).flatten()


    def _jitter_cov_from_params(self, Jx, Jy, Jxy_deg):
        """Convert principal-axis jitter parameters [mas, mas, deg] to covariance."""
        theta = torch.deg2rad(Jxy_deg)
        c, s = torch.cos(theta), torch.sin(theta)
        R = torch.stack((
            torch.stack((c, -s), dim=-1),
            torch.stack((s,  c), dim=-1),
        ), dim=-2)
        D = torch.zeros(*Jx.shape, 2, 2, device=Jx.device, dtype=Jx.dtype)
        D[..., 0, 0] = Jx**2
        D[..., 1, 1] = Jy**2
        return R @ D @ R.transpose(-1, -2)


    def _jitter_params_from_cov_mas(self, Sigma_mas2):
        """Convert covariance in mas² to TipTorch jitter parameters."""
        if cov_to_jitter_mas is not None:
            # cov_to_jitter_mas expects TT covariance in nm² and applies the
            # nm->mas scaling itself, so do not use it for already-mas covariances.
            pass

        eigvals, eigvecs = torch.linalg.eigh(Sigma_mas2)
        order = torch.argsort(eigvals, dim=-1, descending=True)
        eigvals = torch.gather(eigvals, -1, order)
        batch_shape = Sigma_mas2.shape[:-2]
        eigvecs = torch.gather(eigvecs, -1, order[..., None, :].expand(*batch_shape, 2, 2))

        Jx = torch.sqrt(torch.clamp(eigvals[..., 0], min=0.0))
        Jy = torch.sqrt(torch.clamp(eigvals[..., 1], min=0.0))
        Jxy = torch.rad2deg(torch.atan2(eigvecs[..., 1, 0], eigvecs[..., 0, 0]))
        return Jx, Jy, Jxy


    def _static_jitter_params(self):
        """Return telescope static jitter as TipTorch parameters in [mas, mas, deg]."""
        if self.jitter_FWHM is None:
            return None

        sigma_from_fwhm = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))

        if isinstance(self.jitter_FWHM, list):
            Jx = float(self.jitter_FWHM[0]) * sigma_from_fwhm
            Jy = float(self.jitter_FWHM[1]) * sigma_from_fwhm
            # The historical P3 ellipse stored the angle in radians.
            Jxy = np.rad2deg(float(self.jitter_FWHM[2]))
        else:
            Jx = Jy = float(self.jitter_FWHM) * sigma_from_fwhm
            Jxy = 0.0

        n_src = self.nPointings
        dtype = self.model.wvl.dtype
        device = self.model.device
        return (
            torch.full((n_src,), Jx,  device=device, dtype=dtype),
            torch.full((n_src,), Jy,  device=device, dtype=dtype),
            torch.full((n_src,), Jxy, device=device, dtype=dtype),
        )


    def _LO_jitter_params(self):
        """Convert MavisLO tip/tilt covariance to TipTorch jitter parameters."""
        if not (self.LOisOn and self.doConvolve and self.LO_tiptilt_cov is not None):
            return None

        Sigma_tt_nm2 = torch.as_tensor(
            self.LO_tiptilt_cov,
            device=self.model.device,
            dtype=self.model.wvl.dtype,
        )

        if cov_to_jitter_mas is not None:
            return cov_to_jitter_mas(Sigma_tt_nm2, telescope_diameter=2.0*self.tel_radius)

        # Fallback equivalent to utils.cov_to_jitter_mas().
        Jx_nm, Jy_nm, Jxy = self._jitter_params_from_cov_mas(Sigma_tt_nm2)
        scale = (1.0 / rad2mas) * (2.0*self.tel_radius) / (4.0e-9)
        return Jx_nm / scale, Jy_nm / scale, Jxy


    def _combined_jitter_params(self):
        """Combine static telescope jitter and LO residual jitter in covariance space."""
        static = self._static_jitter_params()
        lo = self._LO_jitter_params()

        if static is None and lo is None:
            zeros = torch.zeros(self.nPointings, device=self.model.device, dtype=self.model.wvl.dtype)
            return zeros, zeros, zeros
        if static is None:
            return lo
        if lo is None:
            return static

        if combine_zero_centered_jitters is not None:
            return combine_zero_centered_jitters(*lo, *static)

        # Fallback: convolution of two centered Gaussian jitters => add covariances.
        Sigma = self._jitter_cov_from_params(*lo) + self._jitter_cov_from_params(*static)
        return self._jitter_params_from_cov_mas(Sigma)


    def _compute_reference_tiptorch_psd(self):
        """Compute/cache the HO PSD with TipTorch, without any jitter convolution."""
        if self.verbose:
            print('******** TipTorch HO PSD science directions')

        self._zero_model_jitter()
        with torch.no_grad():
            self.PSD = self.model.ComputePSD().detach().clone()

        self.HO_res = torch.sqrt(torch.clamp(self.PSD.sum(dim=(-2, -1)), min=0.0)).cpu().numpy()
        if self.HO_res.ndim > 1:
            self.HO_res = self.HO_res[:, 0]


    def _source_model(self, zenith, azimuth, wavelength):
        """Create a temporary TipTorch model for guide-star/focus directions."""
        cfg = copy.deepcopy(self.config_torch)
        device = self.model.device
        dtype = self.model.wvl.dtype
        zenith = torch.as_tensor(zenith,  device=device, dtype=dtype).flatten()
        azimuth = torch.as_tensor(azimuth, device=device, dtype=dtype).flatten()
        wavelength = torch.as_tensor([wavelength], device=device, dtype=dtype).flatten()

        cfg['NumberSources'] = int(zenith.numel())
        cfg['sources_science']['Zenith'] = zenith
        cfg['sources_science']['Azimuth'] = azimuth
        cfg['sources_science']['Wavelength'] = wavelength

        model = TipTorch(
            AO_config=cfg,
            norm_regime='sum',
            device=device,
            oversampling=self.overSamp,
            retain_PSDs=True,
        )
        self._zero_model_jitter(model)
        return model


    def _sensor_PSF_sampling(self, sensor_wvl, sensor_pixel_scales):
        """Return the TipTorch PSF sampling for a sensor wavelength."""
        psf_pixel_scale_mas = self.psInMas * sensor_wvl / self.wvlMax
        return psf_pixel_scale_mas, self.nPixPSF, False


    def _compute_sensor_psf_metrics_tiptorch(self, zenith, azimuth, wavelength, pixel_scales, label):
        """Generate guide-star PSFs with TipTorch and measure SR/FWHM/EE."""
        model = self._source_model(zenith, azimuth, wavelength)
        psf_pixel_scale_mas, nPixPSF, _ = self._sensor_PSF_sampling(wavelength, pixel_scales)

        with torch.no_grad():
            PSF = model().detach().cpu().numpy()[:, 0]
            PSD = model.PSD.detach().cpu().numpy()[:, 0]
            PSF_DL = model.DLPSF().detach().cpu().numpy()[0, 0]

        sr_list, fwhm_list, ee_list = [], [], []
        for idx, psf in enumerate(PSF):
            sr = np.max(psf) / np.max(PSF_DL) * PSF_DL.sum() / psf.sum()
            fwhm = getFWHM(psf, psf_pixel_scale_mas, method='contour', nargout=1)

            sensor_pixel_scale_i = pixel_scales[idx] if isinstance(pixel_scales, list) else pixel_scales
            if 2 * fwhm >= nPixPSF * psf_pixel_scale_mas:
                ee = 1.0
            else:
                ee_curve, rr = getEncircledEnergy(
                    psf,
                    pixelscale=psf_pixel_scale_mas,
                    center=centeredPixelCoords(nPixPSF),
                    nargout=2,
                )
                ee_curve *= 1 / np.max(ee_curve)
                ee = interp1d(rr, ee_curve, kind='cubic', bounds_error=False)(max(fwhm, sensor_pixel_scale_i))

            sr_list.append(float(sr))
            fwhm_list.append(float(fwhm))
            ee_list.append(float(np.asarray(ee)))

            if self.verbose:
                print('SR(@', int(wavelength * 1e9), 'nm)        :', "%.5f" % sr)
                print('FWHM(@', int(wavelength * 1e9), 'nm) [mas]:', "%.3f" % fwhm)
                print(label, ':', "%.5f" % float(np.asarray(ee)))

        return sr_list, fwhm_list, ee_list


    def ngsPSF(self):
        """Compute guide-star PSF quality inputs for MavisLO using TipTorch."""
        if self.verbose:
            print('******** TipTorch LO PSF - NGS directions')

        self.NGS_SR_field, self.NGS_FWHM_mas_field, self.NGS_EE_field = \
            self._compute_sensor_psf_metrics_tiptorch(
                self.LO_zen_field,
                self.LO_az_field,
                self.LO_wvl,
                self.LO_psInMas,
                label='EE                  '
            )

        if self.addLOAlias:
            if self.verbose:
                print('Adding aliasing error on LO!')
            dl_fwhm = []
            model = self._source_model(self.LO_zen_field, self.LO_az_field, self.LO_wvl)
            with torch.no_grad():
                dl = model.DLPSF().detach().cpu().numpy()
            psf_pixel_scale_mas, _, _ = self._sensor_PSF_sampling(self.LO_wvl, self.LO_psInMas)
            for psf in dl[:, 0]:
                dl_fwhm.append(float(getFWHM(psf, psf_pixel_scale_mas, method='contour', nargout=1)))
            self.NGS_DL_FWHM_mas = dl_fwhm
        else:
            self.NGS_DL_FWHM_mas = None

        if not self.addFocusError:
            return

        if 'sensor_Focus' not in self.config:
            self.Focus_SR_field = self.NGS_SR_field
            self.Focus_FWHM_mas_field = self.NGS_FWHM_mas_field
            self.Focus_EE_field = self.NGS_EE_field
            return

        if self.verbose:
            print('******** TipTorch Focus Sensor PSF - NGS directions')

        self.Focus_SR_field, self.Focus_FWHM_mas_field, self.Focus_EE_field = \
            self._compute_sensor_psf_metrics_tiptorch(
                self.LO_zen_field,
                self.LO_az_field,
                self.Focus_wvl,
                self.Focus_psInMas,
                label='EE (focus sensor)    '
            )


    def _prepare_static_PSF_state(self, astIndex):
        """Prepare TipTorch HO state and MavisLO guide-star inputs."""
        if not (astIndex is None or self.firstSimCall):
            return

        self._compute_reference_tiptorch_psd()

        if not self.LOisOn:
            return

        if self.verbose:
            print('******** LO PART')

        self.ngsPSF()
        self.mLO = MavisLO(self.path, self.parametersFile, verbose=self.verbose)


    def _focus_filter_full_psd(self):
        """Build a normalized full-plane TipTorch focus filter in PSD coordinates."""
        focus_half = self.model.FocusFilter(self.model.k2).real
        focus_full = self.model.half_PSD_to_full(focus_half).real
        return focus_full / torch.clamp(focus_full.sum(dim=(-2, -1), keepdim=True), min=1e-30)


    def _compute_full_field_LO_terms(self):
        """Compute full-field MavisLO tip/tilt and optional focus terms."""
        self.LO_tiptilt_cov = self.mLO.computeTotalResidualMatrix(
            np.array(self.cartSciencePointingCoords),
            self.cartNGSCoords_field,
            self.NGS_fluxes_field,
            self.LO_freqs_field,
            self.NGS_SR_field,
            self.NGS_EE_field,
            self.NGS_FWHM_mas_field,
            aNGS_FWHM_DL_mas=self.NGS_DL_FWHM_mas,
            doAll=True,
        )
        self.LO_res = np.sqrt(np.trace(self.LO_tiptilt_cov, axis1=1, axis2=2))

        if not self.addFocusError:
            return

        self.focus_cov = self.mLO.computeFocusTotalResidualMatrix(
            self.cartNGSCoords_field,
            self.Focus_fluxes_field,
            self.Focus_freqs_field,
            self.Focus_SR_field,
            self.Focus_EE_field,
            self.Focus_FWHM_mas_field,
        )
        self.GF_res = np.sqrt(max(self.focus_cov[0], 0))
        self.GFinPSD = True

        # Preserve the old full-field convention: focus is folded into the HO PSD.
        self.PSD = self.PSD + (self.GF_res**2) * self._focus_filter_full_psd()


    def _prepare_asterism_LO_cache(self):
        """Initialize MavisLO caches before indexed per-asterism calls."""
        if not self.firstSimCall:
            return

        self.mLO.computeTotalResidualMatrix(
            np.array(self.cartSciencePointingCoords),
            self.cartNGSCoords_field,
            self.NGS_fluxes_field,
            self.LO_freqs_field,
            self.NGS_SR_field,
            self.NGS_EE_field,
            self.NGS_FWHM_mas_field,
            aNGS_FWHM_DL_mas=self.NGS_DL_FWHM_mas,
            doAll=False,
        )

        if self.addFocusError:
            self.mLO.computeFocusTotalResidualMatrix(
                self.cartNGSCoords_field,
                self.Focus_fluxes_field,
                self.Focus_freqs_field,
                self.Focus_SR_field,
                self.Focus_EE_field,
                self.Focus_FWHM_mas_field,
            )


    def finalPSF(self, astIndex=None):
        """Generate final science PSFs with TipTorch.

        The final jitter is not applied by an external image convolution anymore.
        Instead, LO residual covariance and static telescope jitter are converted
        to a single Gaussian jitter ellipse and passed to TipTorch through Jx/Jy/Jxy.
        TipTorch then applies it as an OTF-domain convolution kernel.
        """
        if self.verbose:
            print('******** TipTorch FINAL PSF')

        if self.LOisOn and self.LO_tiptilt_cov is not None:
            self.cov_ellipses = np.column_stack([
                self._LO_jitter_params()[2].detach().cpu().numpy(),
                self._LO_jitter_params()[0].detach().cpu().numpy(),
                self._LO_jitter_params()[1].detach().cpu().numpy(),
            ])
            if not self.doConvolveAsterism:
                return

        Jx, Jy, Jxy = self._combined_jitter_params() if self.doConvolve else (
            torch.zeros(self.nPointings, device=self.model.device, dtype=self.model.wvl.dtype),
            torch.zeros(self.nPointings, device=self.model.device, dtype=self.model.wvl.dtype),
            torch.zeros(self.nPointings, device=self.model.device, dtype=self.model.wvl.dtype),
        )
        self._set_model_jitter(Jx, Jy, Jxy)

        with torch.no_grad():
            self.PSF_tensor = self.model(PSD=self.PSD).detach()

        self.cubeResultsArray = self.PSF_tensor.detach().cpu().numpy()
        self.results = self.cubeResultsArray
        self.cubeResults = self.cubeResultsArray

        self.pointings_FWHM_mas = []
        for i, wvl in enumerate(self.wvl):
            fwhm_list = []
            for psf in self.cubeResultsArray[:, i]:
                fwhm_list.append(float(getFWHM(psf, self.psInMas, method='contour', nargout=1)))
            self.pointings_FWHM_mas.append(fwhm_list)
        if self.nWvl == 1:
            self.pointings_FWHM_mas = self.pointings_FWHM_mas[0]

        self.final_Jx_mas = Jx.detach().cpu().numpy()
        self.final_Jy_mas = Jy.detach().cpu().numpy()
        self.final_Jxy_deg = Jxy.detach().cpu().numpy()


    def _plot_final_PSFs(self):
        if not self.doPlot:
            return
        tiledDisplay(self.cubeResultsArray[:, 0])


    def _finalize_full_field_results(self, astIndex):
        """Compute TipTorch open-loop and diffraction-limited reference PSFs."""
        if astIndex is not None:
            return
        with torch.no_grad():
            self.psfOL = self.model.OLPSF(include_static=True, include_jitter=False).detach().cpu().numpy()
            self.psfDL = self.model.DLPSF().detach().cpu().numpy()

        if self.verbose:
            print('HO_res [nm]:', self.HO_res)
            if self.LOisOn and hasattr(self, 'LO_res'):
                print('LO_res [nm]:', self.LO_res)
            if hasattr(self, 'GF_res'):
                print('GF_res [nm]:', self.GF_res)


    def doOverallSimulation(self, astIndex=None):
        if self.LOisOn:
            self.configLO()

        self.results = []

        self._prepare_static_PSF_state(astIndex)
        self._compute_LO_terms(astIndex)
        self.finalPSF(astIndex)
        self._plot_final_PSFs()

        self.firstSimCall = False
        self._finalize_full_field_results(astIndex)


    '''
    def computePSF1D(self):
        psf1d = []
        psf1d_radius = None
        psf1d_radius_list_list = []
        # === Precompute polar grid once if SupSamp flag ===
        use_polar_interp = self.SupSamp and self.SupSamp[1] == 2
        polar_grid = None
        r_vals_interp = None
        
        if use_polar_interp:
            step_interp = self.SupSamp[0]
            first_psf = self.cubeResults[0][0] if self.nWvl > 1 else self.cubeResults[0]
            center = np.unravel_index(np.argmax(first_psf), first_psf.shape)
            maxradius = self.psInMas * (first_psf.shape[0] / 2)
            r_vals_interp, polar_grid = precompute_polar_grid(step=step_interp, pixelscale=self.psInMas,
                                                              maxradius=maxradius, center=center)
        for i in range(self.nWvl):
            if self.nWvl>1:
                cubeResults = self.cubeResults[i]
            else:
                cubeResults = self.cubeResults
            psf1dList= []
            psf1d_radius_list = []
            for psf in cubeResults:
                psfRadius = psf.shape[0]/2
                center = np.unravel_index(np.argmax(psf), psf.shape)
                rr, radialprofile, ee = radial_profile(psf,
                                                       ext=0,
                                                       pixelscale=self.psInMas,
                                                       ee=True,
                                                       center=center,
                                                       stddev=False,
                                                       binsize=None,
                                                       maxradius=self.psInMas*psfRadius,
                                                       normalize='total',
                                                       pa_range=None,
                                                       slice=0,
                                                       nargout=2,
                                                       supersamp=self.SupSamp, 
                                                       polar_grid=polar_grid,
                                                       r_vals=r_vals_interp,
                                                       verbose=self.verbose)
                psf1dList.append(radialprofile)
                psf1d_radius = rr
                psf1d_radius_list.append(rr)
            psf1d.append(psf1dList)
            psf1d_radius_list_list.append(psf1d_radius_list)
        self.psf1d = np.asarray(psf1d)
        self.psf1d_radius = np.asarray(psf1d_radius)
        self.psf1d_radius_list_list = np.asarray(psf1d_radius_list_list)
        self.psf1d_data = np.vstack( (self.psf1d_radius_list_list, self.psf1d) )


    def savePSFprofileJSON(self):
        now = datetime.now()
        psf_data = {}
        psf_data['radius'] = self.psf1d_radius.tolist()
        psf_data['psf'] = self.psf1d.tolist()
        filename = os.path.join(self.outputDir, self.outputFile + '1D_PSF' + '.json')
        execution_infos = {}
        execution_infos['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        execution_infos['TIPTOP version'] = __version__
        jsondict = {}
        jsondict['execution_infos'] = execution_infos
        jsondict['infos'] = self.config
        jsondict['psf'] = psf_data
        with open(filename, 'w') as f:
            json.dump(jsondict, f)


    def saveResults(self):
        # save PSF cube in fits
        hdul1 = fits.HDUList()
        hdul1.append(fits.PrimaryHDU())
        hdul1.append(fits.ImageHDU(data=self.cubeResultsArray))
        hdul1.append(fits.ImageHDU(data=cpuArray(self.psfOL.sampling))) # append open-loop PSF
        hdul1.append(fits.ImageHDU(data=cpuArray(self.psfDL.sampling))) # append diffraction limited PSF
        if self.savePSDs:
            hdul1.append(fits.ImageHDU(data=cpuArray(self.PSD))) # append high order PSD
        hdul1.append(fits.ImageHDU(data=cpuArray(self.psf1d_data))) # append radial profiles forthe final PSFs

        now = datetime.now()
        # header
        hdr0 = hdul1[0].header       
        hdr0['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        hdr0['TIPTOP_V'] = __version__
        # parameters in the header
        for key_primary in self.config:
            for key_secondary in self.config[key_primary]:
                temp = self.config[key_primary][key_secondary]
                if isinstance(temp, list):
                    iii = 0
                    for elem in temp:
                        if isinstance(elem, list):
                            jjj = 0
                            for elem2 in elem:
                                add_hdr_keyword(hdr0,key_primary,key_secondary,elem2,iii=str(iii),jjj=str(jjj))
                                jjj += 1
                        else:                        
                            add_hdr_keyword(hdr0,key_primary,key_secondary,elem,iii=str(iii))
                        iii += 1
                else:
                    add_hdr_keyword(hdr0, key_primary,key_secondary,temp)

        # header of the PSFs
        hdr1 = hdul1[1].header
        hdr1['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        hdr1['CONTENT'] = "PSF CUBE"
        hdr1['SIZE'] = str(self.cubeResultsArray.shape)
        if self.nWvl>1:
            for i in range(self.nWvl):
                hdr1['WL_NM'+str(i).zfill(3)] = str(int(self.wvl[i]*1e9))
        else:
            hdr1['WL_NM'] = str(int(self.wvl[0]*1e9))
        hdr1['PIX_MAS'] = str(self.psInMas)
        hdr1['CC'] = "CARTESIAN COORD. IN ASEC OF THE "+str(self.pointings.shape[1])+" SOURCES"
        
        for i in range(self.pointings.shape[1]):
            hdr1['CCX' + str(i).zfill(4)] = np.round(self.pointings[0, i], 3).item()
            hdr1['CCY' + str(i).zfill(4)] = np.round(self.pointings[1, i], 3).item()
            
        if hasattr(self,'HO_res'):
            hdr1['RESH'] = "High Order residual in nm RMS"
            for i in range(self.HO_res.shape[0]):
                hdr1['RESH'+str(i).zfill(4)] =  np.round(cpuArray(self.HO_res[i]), 3)
                
        if hasattr(self,'LO_res'):
            hdr1['RESL'] = "Low Order residual in nm RMS"
            for i in range(self.LO_res.shape[0]):
                hdr1['RESL'+str(i).zfill(4)] = np.round(cpuArray(self.LO_res[i]), 3)
                
        if hasattr(self,'GF_res'):
            hdr1['RESF'] = "Global Focus residual in nm RMS (included in PSD)"
            hdr1['RESF0000'] = np.round(cpuArray(self.GF_res), 3)
            
        if self.addSrAndFwhm:
            for i in range(self.nWvl):
                if self.nWvl>1:
                    cubeResultsArray = self.cubeResultsArray[i]
                    wTxt = 'W'+str(i).zfill(2)
                    fTxt = 'FW'
                    eTxt = 'EE'
                    Nfill = 2
                else:
                    cubeResultsArray = self.cubeResultsArray
                    wTxt = ''
                    fTxt = 'FWHM'
                    eTxt = 'EE'+str(np.round(self.eeRadiusInMas))
                    Nfill = 4
                samp = self.wvl[i] * rad2mas / (self.psInMas*2*self.tel_radius)
                
                for j in range(cubeResultsArray.shape[0]):
                    sr_temp = getStrehl(cubeResultsArray[j,:,:], self.fao.ao.tel.pupil,
                                        samp, method='max', psfInOnePix=True)
                    hdr1['SR'+str(j).zfill(Nfill)+wTxt] = float(np.round(sr_temp,5))
                    
                for j in range(cubeResultsArray.shape[0]):
                    fwhm_temp = getFWHM(cubeResultsArray[j,:,:], self.psInMas, method='contour', nargout=1)
                    hdr1[fTxt+str(j).zfill(Nfill)+wTxt] = np.round(fwhm_temp, 3)
                    
                for j in range(cubeResultsArray.shape[0]):
                    if self.ensquaredEnergy:
                        ee = cpuArray(getEnsquaredEnergy(cubeResultsArray[j,:,:]))
                        rr = np.arange(1, ee.shape[0]*2, 2) * self.psInMas * 0.5
                    else:
                        ee,rr = getEncircledEnergy(cubeResultsArray[j,:,:], pixelscale=self.psInMas, center=centeredPixelCoords(self.nPixPSF), nargout=2)
                    ee_at_radius_fn = interp1d(rr, ee, kind='cubic', bounds_error=False)
                    hdr1[eTxt+str(j).zfill(Nfill)+wTxt] = np.round(ee_at_radius_fn(self.eeRadiusInMas).take(0),5)

        # header of the OPEN-LOOP PSF
        hdr2 = hdul1[2].header
        hdr2['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        hdr2['CONTENT'] = "OPEN-LOOP PSF"
        hdr2['SIZE'] = str(self.psfOL.sampling.shape)

        # header of the DIFFRACTION LIMITED PSF
        hdr3 = hdul1[3].header
        hdr3['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        hdr3['CONTENT'] = "DIFFRACTION LIMITED PSF"
        hdr3['SIZE'] = str(self.psfDL.sampling.shape)

        ii = 4
        if self.savePSDs:
            # header of the PSD
            hdr4 = hdul1[4].header
            hdr4['TIME'] = now.strftime("%Y%m%d_%H%M%S")
            hdr4['CONTENT'] = "High Order PSD"
            hdr4['SIZE'] = str(self.PSD.shape)
            ii = 5

        # header of the Total PSFs profiles
        hdr5 = hdul1[ii].header
        hdr5['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        hdr5['CONTENT'] = "Final PSFs profiles"
        hdr5['SIZE'] = str(self.psf1d_data.shape)
        if self.SupSamp:
            hdr5['SAMP_MAS'] = str(self.SupSamp[0])

        hdul1.writeto( os.path.join(self.outputDir, self.outputFile + '.fits'), overwrite=True)
        if self.verbose:
            print("Output cube shape:", self.cubeResultsArray.shape)
            print("Output dtype:", self.cubeResultsArray.dtype)


    def computeOL_PSD(self):
        # OPEN-LOOP PSD
        k = np.sqrt(self.fao.freq.k2_)
        pf = FourierUtils.pistonFilter(2*self.tel_radius,k)
        spectrum = arrayP3toMastsel(self.fao.ao.atm.spectrum(k) * pf)
        psdOL = Field(self.wvlRef, self.N, self.freq_range, 'rad')
        psdOL.sampling = spectrum * (self.dk*self.wvlRef/np.pi)**2 # the PSD must be provided in m^2.m^2
        padPSD = self.nWvl > 1
        mask  = arrayP3toMastsel(self.fao.ao.tel.pupil)
        psfOL = psdSetToPsfSet([psdOL.sampling], mask,
                                self.wvlRef, self.N, self.sx, self.grid_diameter,
                                self.freq_range, self.dk, self.nPixPSF,
                                self.wvlMax, self.overSamp, padPSD=padPSD)
        self.psfOL = psfOL[0]
        
        if self.doPlot:
            fig, ax1 = plt.subplots(1,1)
            im = ax1.imshow(np.log(np.abs(cpuArray(self.psfOL.sampling)) + 1e-20), cmap='hot')
            ax1.set_title('open loop PSF', color='black')


    def computeDL_PSD(self):
        # DIFFRACTION LIMITED PSD
        psdDL = Field(self.wvlRef, self.N, self.freq_range, 'rad')
        padPSD = self.nWvl > 1
        mask = arrayP3toMastsel(self.fao.ao.tel.pupil)
        psfDL = psdSetToPsfSet([psdDL.sampling], mask,
                                self.wvlRef, self.N, self.sx, self.grid_diameter,
                                self.freq_range, self.dk, self.nPixPSF,
                                self.wvlMax, self.overSamp, padPSD=padPSD)
        self.psfDL = psfDL[0]
        if self.doPlot:
            fig, ax2 = plt.subplots(1,1)
            im = ax2.imshow(np.log(np.abs(cpuArray(self.psfDL.sampling)) + 1e-20), cmap='hot')
            ax2.set_title('diffraction limited PSF', color='black')
'''
