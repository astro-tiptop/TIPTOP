import numpy as np
from matplotlib import rc
from scipy.interpolate import interp1d

# P3 Imports
from p3.aoSystem.fourierModel import fourierModel
from p3.aoSystem.FourierUtils import getStrehl, getFWHM, getEncircledEnergy, getEnsquaredEnergy, pistonFilter

# Mastsel Imports
from mastsel import MavisLO, psdSetToPsfSet, longExposurePsf, Field, convolve, residualToSpectrum, maskSA
from mastsel.mavisPsf import centeredPixelCoords, padOrCropCentered, mastselPsfPrecision
from mastsel.mavisUtilities import sigma_from_FWHM, tiledDisplay, plotEllipses, polarToCartesian, congrid

# Tiptop Imports
from .tiptopUtils import arrayP3toMastsel, cpuArray
from .abstractSimulation import AbstractSimulation

rc("text", usetex=False)
rad2mas = 3600 * 180 * 1000 / np.pi

class baseSimulation(AbstractSimulation):
    """
    P3-based Backend Implementation for TIPTOP Simulations.
    Inherits the strict lifecycle template from AbstractSimulation.
    """

    def __init__(self, path, parametersFile, outputDir, outputFile, doConvolve=True,
                 doPlot=False, addSrAndFwhm=True, verbose=False, getHoErrorBreakDown=False,
                 savePSDs=False, ensquaredEnergy=False, eeRadiusInMas=50):
        
        # Superclass handles all configuration loading, validation, and standard properties
        super().__init__(path, parametersFile, outputDir, outputFile, doConvolve,
                         doPlot, addSrAndFwhm, verbose, getHoErrorBreakDown,
                         savePSDs, ensquaredEnergy, eeRadiusInMas)
        
        # P3-specific state variables
        self.fao = None
        self.mLO = None

    # ----------------------------------------------------------------------------
    # --- IMPLEMENTATION OF ABSTRACT METHODS ---
    # ----------------------------------------------------------------------------

    def _configure_LO_parameters(self, astIndex):
        """
        Configures the geometry and low order arrays based on the loaded configuration.
        """
        if not self.LOisOn:
            return

        self.cartSciencePointingCoords = np.dstack(
            (self.pointings[0, :], self.pointings[1, :])
        ).reshape(-1, 2)
        
        # Extract LO wavelength
        LO_wvl_temp = self.my_data_map['sources_LO']['Wavelength']
        self.LO_wvl = LO_wvl_temp[0] if isinstance(LO_wvl_temp, list) else LO_wvl_temp

        self.LO_zen_field = self.my_data_map['sources_LO']['Zenith']
        self.LO_az_field = self.my_data_map['sources_LO']['Azimuth']
        self.LO_fluxes_field = self.my_data_map['sensor_LO']['NumberPhotons']
        
        # Format Pixel Scale
        lo_ps = self.my_data_map['sensor_LO']['PixelScale']
        self.LO_psInMas = lo_ps if isinstance(lo_ps, list) else [lo_ps] * len(self.LO_zen_field)
        
        # Format Frame Rate
        lo_fr = self.my_data_map['RTC']['SensorFrameRate_LO']
        self.LO_freqs_field = lo_fr if isinstance(lo_fr, list) else [lo_fr] * len(self.LO_zen_field)
        
        self.addLoAlias = self.my_data_map['sensor_LO'].get('addAliasError', False)

        # Focus Sensor Configuration
        if 'sensor_Focus' in self.my_data_map:
            self.Focus_fluxes4s_field = self.my_data_map['sensor_Focus']['NumberPhotons']
            f_ps = self.my_data_map['sensor_Focus']['PixelScale']
            self.Focus_psInMas = f_ps if isinstance(f_ps, list) else [f_ps] * len(self.LO_zen_field)
            
            f_wvl = self.my_data_map.get('sources_Focus', {}).get('Wavelength', self.my_data_map['sources_LO']['Wavelength'])
            self.Focus_wvl = f_wvl[0] if isinstance(f_wvl, list) else f_wvl
        else:
            self.Focus_fluxes4s_field = self.LO_fluxes_field
            self.Focus_psInMas = self.LO_psInMas
            self.Focus_wvl = self.LO_wvl

        f_reqs = self.my_data_map['RTC'].get('SensorFrameRate_Focus', self.LO_freqs_field)
        self.Focus_freqs_field = f_reqs if isinstance(f_reqs, list) else [f_reqs] * len(self.LO_zen_field)

        self.NGS_fluxes_field = [f * fr for f, fr in zip(self.LO_fluxes_field, self.LO_freqs_field)]
        self.Focus_fluxes_field = [f * fr for f, fr in zip(self.Focus_fluxes4s_field, self.Focus_freqs_field)]
        
        # Coordinate arrays
        polarNGSCoords = np.column_stack((self.LO_zen_field, self.LO_az_field))
        self.nNaturalGS_field = len(self.LO_zen_field)
        self.cartNGSCoords_field = np.asarray([polarToCartesian(polarNGSCoords[i, :]) for i in range(self.nNaturalGS_field)])
        self.currentAsterismIndices = list(range(len(self.LO_zen_field)))
        
        self._set_asterism_data()

    def _set_asterism_data(self):
        """Helper to set current asterism properties based on indices."""
        self.LO_zen_asterism = [self.LO_zen_field[i] for i in self.currentAsterismIndices]
        self.LO_az_asterism = [self.LO_az_field[i] for i in self.currentAsterismIndices]
        self.LO_fluxes_asterism = [self.LO_fluxes_field[i] for i in self.currentAsterismIndices]
        self.LO_freqs_asterism = [self.LO_freqs_field[i] for i in self.currentAsterismIndices]
        self.NGS_fluxes_asterism = [self.NGS_fluxes_field[i] for i in self.currentAsterismIndices]
        self.Focus_fluxes_asterism = [self.Focus_fluxes_field[i] for i in self.currentAsterismIndices]
        self.cartNGSCoords_asterism = [self.cartNGSCoords_field[i] for i in self.currentAsterismIndices]

    def _prepare_static_PSF_state(self, astIndex):
        """
        Initializes the Fourier HO model and computes the static OL/DL representations.
        """
        if not (astIndex is None or self.firstSimCall):
            return

        if self.verbose:
            print('******** HO PSD science and NGSs directions')

        # Instantiate P3 model
        self.fao = fourierModel(self.fullPathFilename, calcPSF=False, verbose=self.verbose,
                                display=False, getPSDatNGSpositions=self.LOisOn,
                                computeFocalAnisoCov=False, TiltFilter=self.LOisOn,
                                getErrorBreakDown=self.getHoErrorBreakDown, doComputations=False,
                                psdExpansion=True, reduce_memory=True)

        if 'sensor_LO' in self.my_data_map:
            self.fao.my_data_map['sensor_LO']['NumberPhotons'] = self.my_data_map['sensor_LO']['NumberPhotons']
            self.fao.ao.my_data_map['sensor_LO']['NumberPhotons'] = self.my_data_map['sensor_LO']['NumberPhotons']
        if 'sources_LO' in self.my_data_map:
            self.fao.my_data_map['sources_LO'] = self.my_data_map['sources_LO']
            self.fao.ao.my_data_map['sources_LO'] = self.my_data_map['sources_LO']
            self.fao.ao.configLOsensor()
            self.fao.ao.configLO()
            self.fao.ao.configLO_SC()

        if self.verbose:
            print('Setting MASTSEL PSF precision to:', self.fao.dtype)

        mastselPsfPrecision(dtype=self.fao.dtype)

        self.fao.initComputations()

        # Cache geometry needed downstream
        self.PSD = self.fao.PSD.transpose()
        self.N = self.PSD[0].shape[0]
        self.nPointings = self.pointings.shape[1]
        self.nPixPSF = int(self.fao.ao.cam.fovInPix)
        self.overSamp = getattr(self.fao.freq, 'kRef_float', int(self.fao.freq.kRef_))
        self.PSDstep = self.fao.freq.PSDstep
        self.overSamp_lo = self.fao.freq.kGrid_
        self.freq_range = self.N * self.PSDstep
        self.grid_diameter = 1 / self.PSDstep
        self.sx = self.fao.freq.kGrid_
        self.dk = self.fao.freq.dk_
        self.wvlRef = self.fao.freq.wvlRef

        # Setup Mask
        self.mask = Field(self.wvlRef, self.N, self.grid_diameter)
        self.mask.sampling = congrid(arrayP3toMastsel(self.fao.ao.tel.pupil), [self.sx, self.sx])
        self.mask.sampling = padOrCropCentered(self.mask.sampling, self.N, xp=self.mask.xp)

        if abs(float(self.psInMas) - float(cpuArray(self.fao.freq.psInMas[0]))) > 1e-6:
            raise ValueError(f"sensor_science.PixelScale '{self.psInMas}'"
                             f" differs from P3 '{float(cpuArray(self.fao.freq.psInMas[0]))}'")

        self.opdMap = arrayP3toMastsel(self.fao.ao.tel.opdMap_on) if self.fao.ao.tel.opdMap_on is not None else None

        if self.verbose:
            print('PSD step:', self.PSDstep)
            print('PSD freq range:', self.freq_range)
            print('oversampling:', self.overSamp)
            print('sensor_science.PixelScale:', self.psInMas)

        if self.LOisOn:
            if self.verbose:
                print('******** LO PART')
            self._compute_ngs_psf()
            self.mLO = MavisLO(verbose=self.verbose, config_dict=self.my_data_map)

    def _compute_LO_terms(self, astIndex) -> dict:
        """
        PURE FUNCTION: Returns a dictionary of computed terms instead of mutating state.
        """
        if not self.LOisOn:
            return {}

        lo_data = {}

        if astIndex is None:
            # Full Field Processing
            Ctot = self.mLO.computeTotalResidualMatrix(
                np.array(self.cartSciencePointingCoords), self.cartNGSCoords_field, 
                self.NGS_fluxes_field, self.LO_freqs_field, self.NGS_SR_field, 
                self.NGS_EE_field, self.NGS_FWHM_mas_field,
                aNGS_FWHM_DL_mas=self.NGS_DL_FWHM_mas, doAll=True)

            if self.addFocusError:
                CtotFocus = self.mLO.computeFocusTotalResidualMatrix(
                    self.cartNGSCoords_field, self.Focus_fluxes_field,
                    self.Focus_freqs_field, self.Focus_SR_field,
                    self.Focus_EE_field, self.Focus_FWHM_mas_field)
                
                lo_data['GF_res'] = float(np.sqrt(np.maximum(cpuArray(CtotFocus).ravel()[0], 0.0)))
                
                # Apply Global Focus filtering to PSD
                FocusFilter = self.fao.FocusFilter()
                FocusFilter *= 1 / FocusFilter.sum()
                for PSDho in self.PSD:
                    PSDho += (lo_data['GF_res']**2) * FocusFilter
                lo_data['GFinPSD'] = True
        else:
            # Asterism Specific Processing
            if self.firstSimCall:
                self.mLO.computeTotalResidualMatrix(
                    np.array(self.cartSciencePointingCoords), self.cartNGSCoords_field, 
                    self.NGS_fluxes_field, self.LO_freqs_field, self.NGS_SR_field,
                    self.NGS_EE_field, self.NGS_FWHM_mas_field,
                    aNGS_FWHM_DL_mas=self.NGS_DL_FWHM_mas, doAll=False)
                if self.addFocusError:
                    self.mLO.computeFocusTotalResidualMatrix(
                        self.cartNGSCoords_field, self.Focus_fluxes_field,
                        self.Focus_freqs_field, self.Focus_SR_field,
                        self.Focus_EE_field, self.Focus_FWHM_mas_field)

            # Clean faint guide stars
            if np.min(self.NGS_fluxes_asterism) < 1 and np.max(self.NGS_fluxes_asterism) > 1:
                valid_idx = np.where(np.array(self.NGS_fluxes_asterism) > 1)[0]
                self.NGS_fluxes_asterism = [elem for i, elem in enumerate(self.NGS_fluxes_asterism) if i in valid_idx]
                self.Focus_fluxes_asterism = [elem for i, elem in enumerate(self.Focus_fluxes_asterism) if i in valid_idx]
                self.cartNGSCoords_asterism = [elem for i, elem in enumerate(self.cartNGSCoords_asterism) if i in valid_idx]
                self.currentAsterismIndices = [elem for i, elem in enumerate(self.currentAsterismIndices) if i in valid_idx]

            Ctot = self.mLO.computeTotalResidualMatrixI(
                self.currentAsterismIndices, np.array(self.cartSciencePointingCoords),
                np.array(self.cartNGSCoords_asterism), self.NGS_fluxes_asterism)

            if self.addFocusError:
                CtotFocus = self.mLO.computeFocusTotalResidualMatrixI(
                    self.currentAsterismIndices, np.array(self.cartNGSCoords_asterism),
                    self.Focus_fluxes_asterism)
                lo_data['GF_res'] = float(np.sqrt(np.maximum(cpuArray(CtotFocus).ravel()[0], 0.0)))
                lo_data['GFinPSD'] = False

        # Export computed traces
        lo_data['LO_res'] = np.sqrt(np.trace(Ctot, axis1=1, axis2=2))
        if self.doConvolve:
            lo_data['cov_ellipses'] = self.mLO.ellipsesFromCovMats(Ctot)
        
        # Save reference for internal convolutions
        self._Ctot_temp = Ctot 
        
        return lo_data

    def _generate_final_PSF(self, astIndex, lo_data):
        """
        Convolves PSDs and standardizes the output in self.cubeResultsArray.
        """
        if astIndex is None or self.firstSimCall:
            PSD_HO = arrayP3toMastsel(self.PSD[0:self.nPointings])
            mask = arrayP3toMastsel(self.fao.ao.tel.pupil)

            if self.verbose:
                print('******** HO PSF')
                
            psfLongExpPointingsArr = psdSetToPsfSet(
                inputPSDs=PSD_HO,
                mask=mask,
                wavelength=self.wvl,
                nPixPup=self.sx,
                freq_range=self.freq_range,
                dk=self.dk,
                nPixPsf=self.nPixPSF,
                oversampling=self.overSamp,
                opdMap=self.opdMap)

            # Safely compute HO residuals
            self.HO_res = np.sqrt(np.sum(self.PSD[0:self.nPointings], axis=(1, 2)))

            self.pointings_FWHM_mas = []
            for i in range(self.nWvl):
                psfList = psfLongExpPointingsArr[i] if self.nWvl > 1 else psfLongExpPointingsArr
                wvl_c = self.wvl[i] if self.nWvl > 1 else self.wvl[0]
                # Sampling ratio: PSF pixels per lambda/D at this wavelength
                samp_i = wvl_c * rad2mas / (self.psInMas * 2 * self.tel_radius)
                # When samp < 2 the PSF is undersampled; rebin to Nyquist before measuring
                rebin_i = max(1, int(np.ceil(2.0 / samp_i))) if samp_i < 2.0 else 1
                fwhmList = []
                for idx, img in enumerate(psfList):
                    fwhmX, fwhmY = getFWHM(img.sampling, self.psInMas, method='contour',
                                           rebin=rebin_i, nargout=2)
                    fwhm = np.sqrt(fwhmX * fwhmY)
                    fwhmList.append(fwhm)
                    if self.verbose:
                        s1 = cpuArray(PSD_HO[idx]).sum()
                        sr = np.exp(-s1 * (2*np.pi*1e-9/wvl_c)**2)
                        print(f'SR(@{int(wvl_c*1e9)}nm)        : {sr:.5f}')
                        print(f'FWHM(@{int(wvl_c*1e9)}nm) [mas]: {fwhm:.3f}')
                        
                if self.nWvl > 1:
                    self.pointings_FWHM_mas.append(fwhmList)
                else:
                    self.pointings_FWHM_mas = fwhmList

            self.psfLongExpPointingsArr = psfLongExpPointingsArr

        # ------------------------------------------------------------------------
        # --- Final Convolutions ---
        self.results = []
        if self.LOisOn:
            if self.doConvolve:
                if self.doConvolveAsterism:
                    self._apply_final_convolution(lo_data)
                else:
                    self.cov_ellipses = lo_data.get('cov_ellipses')
            else:
                self._apply_no_convolution()
        else:
            self._apply_jitter_convolution()

        # ------------------------------------------------------------------------
        # FORMAT TO STANDARD ARRAY (Expected by AbstractSimulation)
        cubeResultsArray = []
        for i in range(self.nWvl):
            resList = self.results[i] if self.nWvl > 1 else self.results
            cubeResultsArray.append(np.array([cpuArray(img.sampling) for img in resList]))
            
        if self.nWvl > 1:
            self.cubeResultsArray = np.array(cubeResultsArray)
        else:
            self.cubeResultsArray = np.array(cubeResultsArray[0])

    def _finalize_full_field_results(self, astIndex):
        """
        Creates analytical Open-Loop and Diffraction-Limited PSFs.
        Populates standardized NumPy arrays.
        """
        if astIndex is not None:
            return

        # OPEN-LOOP PSD
        k = np.sqrt(self.fao.freq.k2_)
        pf = pistonFilter(2*self.tel_radius, k)
        spectrum = arrayP3toMastsel(self.fao.ao.atm.spectrum(k) * pf)
        psdOL = Field(self.wvlRef, self.N, self.freq_range, 'rad')
        psdOL.sampling = spectrum * (self.dk*self.wvlRef/np.pi)**2
        
        mask = arrayP3toMastsel(self.fao.ao.tel.pupil)
        
        psfOL = psdSetToPsfSet(
                inputPSDs=[psdOL.sampling],
                mask=mask,
                wavelength=self.wvlRef,
                nPixPup=self.sx,
                freq_range=self.freq_range,
                dk=self.dk,
                nPixPsf=self.nPixPSF,
                oversampling=self.overSamp,
                opdMap=self.opdMap)
        self.psf_ol_array = cpuArray(psfOL[0].sampling)

        # DIFFRACTION LIMITED PSD
        psdDL = Field(self.wvlRef, self.N, self.freq_range, 'rad')
        psfDL = psdSetToPsfSet(
                inputPSDs=[psdDL.sampling],
                mask=mask,
                wavelength=self.wvlRef,
                nPixPup=self.sx,
                freq_range=self.freq_range,
                dk=self.dk,
                nPixPsf=self.nPixPSF,
                oversampling=self.overSamp,
                opdMap=self.opdMap)
        self.psf_dl_array = cpuArray(psfDL[0].sampling)

        if self.savePSDs:
            self.psd_array = cpuArray(self.PSD)

        # MANDATORY: Generate 1D radial profiles
        self.computePSF1D()

        if self.verbose:
            print('HO_res [nm]:', self.HO_res)
            if self.LOisOn:
                print('LO_res [nm]:', self.LO_res)
            if self.GF_res is not None:
                print('GF_res [nm]:', self.GF_res)

    def computeMetrics(self):
        """
        Populates Strehl, FWHM, and EE metrics natively.
        """
        self.penalty, self.sr, self.fwhm, self.ee = [], [], [], []
        
        if len(self.results) == 0:
            # Metric fallback if convolution was skipped
            if self.LOisOn:
                base_pen = np.mean(cpuArray(self.LO_res)**2 + cpuArray(self.HO_res)**2)
                if self.addFocusError and not self.GFinPSD:
                    self.penalty.append(np.sqrt(base_pen + self.GF_res**2))
                else:
                    self.penalty.append(np.sqrt(base_pen))
            else:
                self.penalty.append(np.sqrt(np.mean(cpuArray(self.HO_res)**2)))
                
            self.sr.append(np.exp(-4*np.pi**2 * (self.penalty[-1]**2)/(self.wvlRef*1e9)**2))
            scale = (np.pi/(180*3600*1000) * 2 * self.tel_radius / (4*1e-9))
            fwhms_lo = 2.355 * self.LO_res/scale / np.sqrt(2) if self.LOisOn else 0.0
            
            p_fwhm_mas = np.asarray(self.pointings_FWHM_mas) if self.pointings_FWHM_mas is not None else 0.0
            self.fwhm.append(np.sqrt(fwhms_lo**2 + p_fwhm_mas**2))
            self.ee.append(0)
        else:
            for idx, HO_res in enumerate(cpuArray(self.HO_res)):
                if self.LOisOn:
                    if self.addFocusError and not self.GFinPSD:
                        self.penalty.append(np.sqrt(cpuArray(self.LO_res)[idx]**2 + HO_res**2 + self.GF_res**2))
                    else:
                        self.penalty.append(np.sqrt(cpuArray(self.LO_res)[idx]**2 + HO_res**2))
                else:
                    self.penalty.append(HO_res)

            if self.verbose:
                print(f'EE is computed for a radius of {self.eeRadiusInMas} mas')

            for i in range(self.nWvl):
                results_slice = self.results[i] if self.nWvl > 1 else self.results
                samp = self.wvl[i] * rad2mas / (self.psInMas * 2 * self.tel_radius)
                # When samp < 2 the PSF is undersampled; rebin to Nyquist before measuring FWHM
                rebin_fwhm = max(1, int(np.ceil(2.0 / samp))) if samp < 2.0 else 1

                sr_l, fwhm_l, ee_l = [], [], []
                for img in results_slice:
                    sr_l.append(getStrehl(img.sampling, self.fao.ao.tel.pupil, samp, method='max', psfInOnePix=True))
                    fwhm_l.append(getFWHM(img.sampling, self.psInMas, method='contour',
                                          rebin=rebin_fwhm, nargout=1))
                    
                    if self.ensquaredEnergy:
                        ee_ = cpuArray(getEnsquaredEnergy(img.sampling))
                        rr_ = np.arange(1, ee_.shape[0]*2, 2) * self.psInMas * 0.5
                    else:
                        ee_, rr_ = getEncircledEnergy(img.sampling, pixelscale=self.psInMas,
                                                      center=centeredPixelCoords(self.nPixPSF), nargout=2)
                        
                    ee_fn = interp1d(rr_, ee_, kind='cubic', bounds_error=False)
                    ee_val = ee_fn(self.eeRadiusInMas)
                    ee_l.append(float(ee_val.item() if hasattr(ee_val, 'item') else ee_val))

                if self.nWvl > 1:
                    self.sr.append(sr_l)
                    self.fwhm.append(fwhm_l)
                    self.ee.append(ee_l)
                else:
                    self.sr = sr_l
                    self.fwhm = fwhm_l
                    self.ee = ee_l

    def _plot_final_PSFs(self):
        if not self.doPlot:
            return
        res = self.results[0] if self.nWvl > 1 else self.results
        if self.LOisOn and self.doConvolve:
            tiledDisplay(res)
            plotEllipses(self.cartSciencePointingCoords, self.cov_ellipses, 0.4)
        else:
            res[0].standardPlot(True)

    # ----------------------------------------------------------------------------
    # --- P3-SPECIFIC HELPER METHODS ---
    # ----------------------------------------------------------------------------

    def _compute_ngs_psf(self):
        """Generates analytical PSFs for the NGS directions."""
        LO_PSFsInMas = self.psInMas * self.LO_wvl / self.wvlMax
        
        if LO_PSFsInMas / np.min(self.LO_psInMas) > 1 and self.overSamp_lo > 1:
            LO_PSFsInMas /= self.overSamp_lo
            nPixPSFLO = int(self.overSamp_lo * self.nPixPSF)
            lo_oversampling = 1.0
        else:
            nPixPSFLO = self.nPixPSF
            lo_oversampling = self.overSamp

        k = np.sqrt(self.fao.freq.k2_)

        psdNGS_view = arrayP3toMastsel(self.PSD[-self.nNaturalGS_field:])
        
        psdNGS = []

        nSA = self.my_data_map['sensor_LO']['NumberLenslets']
        maskLO = maskSA(nSA, self.nNaturalGS_field, arrayP3toMastsel(self.fao.ao.tel.pupil))

        for i in range(self.nNaturalGS_field):
            nSAi = nSA[i] if len(nSA) == self.nNaturalGS_field else nSA[0]
            if nSAi != 1:
                pf = pistonFilter(2*self.tel_radius/nSAi, k)
                # The result of the multiplication creates a new tensor on the fly,
                # without modifying psdNGS_view (and thus self.PSD)
                psdNGS.append(psdNGS_view[i] * arrayP3toMastsel(pf))
            else:
                # No change needed, we can append the view directly
                psdNGS.append(psdNGS_view[i])

        if self.verbose:
            print('******** LO PSF - NGS directions (1 sub-aperture)')
            
        psfLE_NGS = psdSetToPsfSet(
                inputPSDs=psdNGS,
                mask=maskLO,
                wavelength=self.LO_wvl,
                nPixPup=self.sx,
                freq_range=self.freq_range,
                dk=self.dk,
                nPixPsf=nPixPSFLO,
                oversampling=lo_oversampling,
                opdMap=self.opdMap)

        self.NGS_SR_field, self.NGS_FWHM_mas_field, self.NGS_EE_field = [], [], []
        for idx, img in enumerate(psfLE_NGS):
            s1 = cpuArray(psdNGS[idx]).sum()
            SR = np.exp(-s1 * (2*np.pi*1e-9/self.LO_wvl)**2)
            self.NGS_SR_field.append(SR)
            
            fwhmX, fwhmY = getFWHM(img.sampling, LO_PSFsInMas, method='contour', nargout=2)
            FWHM = np.sqrt(fwhmX * fwhmY)
            self.NGS_FWHM_mas_field.append(FWHM)
            
            if 2 * FWHM >= nPixPSFLO * LO_PSFsInMas:
                ee_NGS = 1.0
            else:
                ee_, rr_ = getEncircledEnergy(img.sampling, pixelscale=LO_PSFsInMas,
                                              center=centeredPixelCoords(nPixPSFLO), nargout=2)
                ee_ *= 1 / np.max(ee_)
                ee_fn = interp1d(rr_, ee_, kind='cubic', bounds_error=False)
                ps_mas_i = self.LO_psInMas[idx] if isinstance(self.LO_psInMas, list) else self.LO_psInMas
                ee_NGS = ee_fn(max([FWHM, ps_mas_i]))
                ee_NGS = float(ee_NGS.item() if hasattr(ee_NGS, 'item') else ee_NGS)
            self.NGS_EE_field.append(ee_NGS)

        if self.addLoAlias:
            self.NGS_DL_FWHM_mas = [] if isinstance(maskLO, list) else None
            len_nSA = len(nSA) if len(nSA) == self.nNaturalGS_field else 1
            for i in range(len_nSA):
                maskI = maskLO[i] if isinstance(maskLO, list) else maskLO
                psdDL = Field(self.LO_wvl, self.N, self.freq_range, 'rad')
                maskField = Field(self.LO_wvl, self.N, self.grid_diameter)
                maskField.sampling = congrid(maskI, [self.sx, self.sx])
                maskField.sampling = padOrCropCentered(maskField.sampling, self.N, xp=maskField.xp)
                psfNgsDL = longExposurePsf(maskField, psdDL)
                fx, fy = getFWHM(psfNgsDL.sampling, LO_PSFsInMas, method='contour', nargout=2)
                f_val = np.sqrt(fx * fy)
                if self.NGS_DL_FWHM_mas is None:
                    self.NGS_DL_FWHM_mas = f_val
                else:
                    self.NGS_DL_FWHM_mas.append(f_val)
        else:
            self.NGS_DL_FWHM_mas = None

        if self.addFocusError:
            # Handle Focus PSFs analogously
            Focus_PSFsInMas = self.psInMas * self.Focus_wvl / self.wvlMax
            if Focus_PSFsInMas / np.min(self.Focus_psInMas) > 1 and self.overSamp_lo > 1:
                Focus_PSFsInMas /= self.overSamp_lo
                nPixPSFFocus = int(self.overSamp_lo * self.nPixPSF)
                focus_oversampling = 1.0
            else:
                nPixPSFFocus = self.nPixPSF
                focus_oversampling = self.overSamp

            if 'sensor_Focus' in self.my_data_map:
                nSAfocus = self.my_data_map['sensor_Focus']['NumberLenslets']
                psdFocus = arrayP3toMastsel(self.PSD[-self.nNaturalGS_field:])
                maskFocus = maskSA(nSAfocus, self.nNaturalGS_field, arrayP3toMastsel(self.fao.ao.tel.pupil))

                for i in range(self.nNaturalGS_field):
                    nSAfocusI = nSAfocus[i] if len(nSAfocus) == self.nNaturalGS_field else nSAfocus[0]
                    if nSAfocusI != 1:
                        pf = pistonFilter(2*self.tel_radius/nSAfocusI, k)
                        psdFocus[i] = psdFocus[i] * pf

                psfLE_Focus = psdSetToPsfSet(
                                inputPSDs=psdFocus,
                                mask=maskFocus,
                                wavelength=self.Focus_wvl,
                                nPixPup=self.sx,
                                freq_range=self.freq_range,
                                dk=self.dk,
                                nPixPsf=nPixPSFFocus,
                                oversampling=focus_oversampling,
                                opdMap=self.opdMap)

                self.Focus_SR_field, self.Focus_FWHM_mas_field, self.Focus_EE_field = [], [], []
                for idx, img in enumerate(psfLE_Focus):
                    s1 = cpuArray(psdFocus[idx]).sum()
                    self.Focus_SR_field.append(np.exp(-s1 * (2*np.pi*1e-9/self.Focus_wvl)**2))
                    
                    fwhmX, fwhmY = getFWHM(img.sampling, Focus_PSFsInMas, method='contour', nargout=2)
                    FWHM = np.sqrt(fwhmX * fwhmY)
                    self.Focus_FWHM_mas_field.append(FWHM)
                    
                    if 2 * FWHM >= nPixPSFFocus * Focus_PSFsInMas:
                        ee_Focus = 1.0
                    else:
                        ee_, rr_ = getEncircledEnergy(img.sampling, pixelscale=Focus_PSFsInMas,
                                                      center=centeredPixelCoords(nPixPSFFocus), nargout=2)
                        ee_ *= 1 / np.max(ee_)
                        ee_fn = interp1d(rr_, ee_, kind='cubic', bounds_error=False)
                        ps_mas_i = self.Focus_psInMas[idx] if isinstance(self.Focus_psInMas, list) else self.Focus_psInMas
                        ee_Focus = ee_fn(max([FWHM, ps_mas_i]))
                        ee_Focus = float(ee_Focus.item() if hasattr(ee_Focus, 'item') else ee_Focus)
                    self.Focus_EE_field.append(ee_Focus)
            else:
                self.Focus_SR_field = self.NGS_SR_field
                self.Focus_FWHM_mas_field = self.NGS_FWHM_mas_field
                self.Focus_EE_field = self.NGS_EE_field

    def _apply_final_convolution(self, lo_data):
        self.cov_ellipses = lo_data.get('cov_ellipses')
        resSpecList, resSpecListJ = [], []
        
        for ellp in self.cov_ellipses:
            ellp = ellp.astype(self.fao.dtype)
            resSpecList.append(residualToSpectrum(ellp, self.wvlRef, self.nPixPSF, 1/(self.nPixPSF * self.psInMas)))
            
            if self.jitter_FWHM is not None:
                ellpJ = [self.jitter_FWHM[2], sigma_from_FWHM(self.jitter_FWHM[0]), sigma_from_FWHM(self.jitter_FWHM[1])] \
                    if isinstance(self.jitter_FWHM, list) \
                        else [0, sigma_from_FWHM(self.jitter_FWHM), sigma_from_FWHM(self.jitter_FWHM)]
                ellpJ = np.array(ellpJ, dtype=self.fao.dtype)
                resSpecListJ.append(residualToSpectrum(ellpJ, self.wvlRef, self.nPixPSF, 1/(self.nPixPSF * self.psInMas)))
            else:
                resSpecListJ.append(0)

        for i in range(self.nWvl):
            psfList = self.psfLongExpPointingsArr[i] if self.nWvl > 1 else self.psfLongExpPointingsArr
            resultList = []
            for psfLongExp, resSpec, resSpecJ in zip(psfList, resSpecList, resSpecListJ):
                temp = convolve(psfLongExp, resSpec)
                if self.jitter_FWHM is not None:
                    temp = convolve(temp, resSpecJ)
                resultList.append(temp)
                
            if self.nWvl > 1:
                self.results.append(resultList)
            else:
                self.results = resultList

    def _apply_no_convolution(self):
        for i in range(self.nWvl):
            psfList = self.psfLongExpPointingsArr[i] if self.nWvl > 1 else self.psfLongExpPointingsArr
            if self.nWvl > 1:
                self.results.append(list(psfList))
            else:
                self.results = list(psfList)

    def _apply_jitter_convolution(self):
        if self.jitter_FWHM is not None:
            ellpJ = [self.jitter_FWHM[2], sigma_from_FWHM(self.jitter_FWHM[0]), sigma_from_FWHM(self.jitter_FWHM[1])] \
                if isinstance(self.jitter_FWHM, list) \
                    else [0, sigma_from_FWHM(self.jitter_FWHM), sigma_from_FWHM(self.jitter_FWHM)]
            resSpecJ = residualToSpectrum(ellpJ, self.wvlRef, self.nPixPSF, 1/(self.nPixPSF * self.psInMas))
            
        for i in range(self.nWvl):
            psfList = self.psfLongExpPointingsArr[i] if self.nWvl > 1 else self.psfLongExpPointingsArr
            resultList = []
            for psfLongExp in psfList:
                resultList.append(convolve(psfLongExp, resSpecJ) if self.jitter_FWHM is not None else psfLongExp)
                
            if self.nWvl > 1:
                self.results.append(resultList)
            else:
                self.results = resultList
