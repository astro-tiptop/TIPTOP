import os
import ast
import json
import yaml
import numpy as np
from abc import ABC, abstractmethod
from astropy.io import fits
from datetime import datetime

from .tiptopUtils import add_hdr_keyword, cpuArray
from mastsel.mavisUtilities import polarToCartesian
from p3.aoSystem.FourierUtils import precompute_polar_grid, radial_profile
from ._version import __version__

rad2mas = 3600 * 180 * 1000 / np.pi

class AbstractSimulation(ABC):
    """
    Abstract Base Class for TIPTOP Simulations.
    
    This class defines the strict lifecycle of a simulation (Template Method Pattern).
    It manages the shared state (configuration, I/O, metrics computation, and FITS saving)
    while delegating backend-specific physics (P3 or TipTorch) to its subclasses via
    pure-ish abstract methods.
    """

    # ----------------------------------------------------------------------------
    # --- Initialization and Configuration Validation ---
    # ----------------------------------------------------------------------------

    def __init__(self, path, parametersFile, outputDir, outputFile, doConvolve=True,
                 doPlot=False, addSrAndFwhm=True, verbose=False, getHoErrorBreakDown=False,
                 savePSDs=False, ensquaredEnergy=False, eeRadiusInMas=50):
        
        self.firstSimCall = True
        self.verbose = verbose
        if self.verbose: 
            np.set_printoptions(precision=3)
            
        self.doConvolveAsterism = True
        self.pointings_FWHM_mas = None
        
        # State variables
        self.path = path
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

        # Standardized Universal Data Containers (To be populated by subclasses)
        self.HO_res = None
        self.LO_res = None
        self.GF_res = None
        self.cubeResultsArray = None # Must be a numpy array of shape (nWvl, nPointings, N, N) or (nPointings, N, N)
        self.psf_ol_array = None     # Must be a 2D numpy array
        self.psf_dl_array = None     # Must be a 2D numpy array
        self.psd_array = None        # Must be a numpy array (if savePSDs is True)
        self.psf1d_data = None
        
        self.results = []
        self.sr = []
        self.fwhm = []
        self.ee = []
        self.penalty = []

        # Load and validate configuration
        self.loadConfigurationFile()
        self._validate_core_configuration()

    def _parse_config_value(self, value):
        if not isinstance(value, str):
            return value
        value = value.strip()
        try:
            return ast.literal_eval(value)
        except (ValueError, SyntaxError):
            return value

    def loadConfigurationFile(self, path=None, parametersFile=None):
        """Loads configuration from .ini or .yml file and populates self.my_data_map."""
        path = path or self.path
        parametersFile = parametersFile or self.parametersFile

        fullPathFilename_ini = os.path.join(path, parametersFile + '.ini')
        fullPathFilename_yml = os.path.join(path, parametersFile + '.yml')

        if os.path.exists(fullPathFilename_yml):
            self.fullPathFilename = fullPathFilename_yml
            with open(fullPathFilename_yml) as f:
                self.my_data_map = yaml.safe_load(f)
        elif os.path.exists(fullPathFilename_ini):
            import configparser
            self.fullPathFilename = fullPathFilename_ini
            config = configparser.ConfigParser()
            config.optionxform = str
            config.read(fullPathFilename_ini)
            self.my_data_map = {}
            for section in config.sections():
                self.my_data_map[section] = {}
                for name, value in config.items(section):
                    self.my_data_map[section].update({name: self._parse_config_value(value)})
        else:
            raise FileNotFoundError(f'No .yml or .ini ({parametersFile}) can be found in {path}')

    def _validate_core_configuration(self):
        """Validates the standard TIPTOP configuration keys irrespective of the backend."""
        def check_sec(sec):
            if sec not in self.my_data_map:
                raise ValueError(f"The section '{sec}' is missing from the parameter file")
        def check_opt(sec, opt):
            if opt not in self.my_data_map.get(sec, {}):
                raise ValueError(f"'{opt}' is missing from section '{sec}'")

        check_sec('telescope')
        check_opt('telescope', 'TelescopeDiameter')
        self.my_data_map['telescope'].setdefault('glFocusOnNGS', False)

        check_sec('sources_science')
        check_opt('sources_science', 'Wavelength')
        self.my_data_map['sources_science'].setdefault('Zenith', [0.0])
        self.my_data_map['sources_science'].setdefault('Azimuth', [0.0])

        if len(self.my_data_map['sources_science']['Zenith']) != len(self.my_data_map['sources_science']['Azimuth']):
            raise ValueError("'Zenith' and 'Azimuth' in 'sources_science' must have the same length")

        check_sec('sensor_science')
        check_opt('sensor_science', 'PixelScale')

        # Handle sensor_science.Super_Sampling
        if 'Super_Sampling' not in self.my_data_map['sensor_science']:
            self.my_data_map['sensor_science']['Super_Sampling'] = None
        else:
            SupSamp_val = self.my_data_map['sensor_science']['Super_Sampling']
            if isinstance(SupSamp_val, (int, float)):
                self.my_data_map['sensor_science']['Super_Sampling'] = [float(SupSamp_val), 2]
            elif isinstance(SupSamp_val, (list, tuple)) and len(SupSamp_val) == 1:
                self.my_data_map['sensor_science']['Super_Sampling'] = [float(SupSamp_val[0]), 2]
            elif isinstance(SupSamp_val, (list, tuple)) and len(SupSamp_val) == 2:
                if int(SupSamp_val[1]) not in (1, 2):
                    raise ValueError("Second value of Super_Sampling must be 1 (1D) or 2 (2D).")
            else:
                raise ValueError("Super_Sampling must be a scalar or list of one/two values.")

        self.LOisOn = 'sensor_LO' in self.my_data_map
        if self.LOisOn:
            check_sec('sources_LO')
            check_opt('sources_LO', 'Wavelength')
            self.my_data_map['sources_LO'].setdefault('Zenith', [0.0])
            self.my_data_map['sources_LO'].setdefault('Azimuth', [0.0])
            check_opt('sensor_LO', 'NumberPhotons')
            check_sec('RTC')
            check_opt('RTC', 'SensorFrameRate_LO')

        # Global Properties extraction
        self.tel_radius = self.my_data_map['telescope']['TelescopeDiameter'] / 2.0
        wvl_temp = self.my_data_map['sources_science']['Wavelength']
        
        if isinstance(wvl_temp, list):
            self.wvlMax = max(wvl_temp)
            self.wvl = wvl_temp
            self.nWvl = len(wvl_temp)
        else:
            self.wvlMax = wvl_temp
            self.wvl = [wvl_temp]
            self.nWvl = 1
            
        self.zenithSrc = self.my_data_map['sources_science']['Zenith']
        self.azimuthSrc = self.my_data_map['sources_science']['Azimuth']
        self.pointings = polarToCartesian(np.array([self.zenithSrc, self.azimuthSrc]))
        self.psInMas = self.my_data_map['sensor_science']['PixelScale']
        self.SupSamp = self.my_data_map['sensor_science']['Super_Sampling']
        
        self.jitter_FWHM = self.my_data_map['telescope'].get('jitter_FWHM', None)
        self.addFocusError = self.my_data_map['telescope'].get('glFocusOnNGS', False)
        self.GFinPSD = False

    # ----------------------------------------------------------------------------
    # --- TEMPLATE METHOD (The Core Architecture) ---
    # ----------------------------------------------------------------------------

    def doOverallSimulation(self, astIndex=None):
        """
        TEMPLATE METHOD: Orchestrates the entire simulation lifecycle.
        It calls abstract methods that MUST be implemented by subclasses (P3 or TipTorch).
        This guarantees strict architectural consistency across backends.
        """
        # 1. Config LO parameters (if active)
        self._configure_LO_parameters(astIndex)
        
        # 2. Prepare the static/High-Order backend state
        self._prepare_static_PSF_state(astIndex)
        
        # 3. Compute Low-Order terms purely (in -> out), preventing temporal coupling
        lo_data = self._compute_LO_terms(astIndex)
        
        # Unpack the returned dictionary to maintain backward compatibility with external wrappers
        self._unpack_LO_data(lo_data)
        
        # 4. Generate Final PSFs
        self._generate_final_PSF(astIndex, lo_data)
        
        # 5. Optional Plotting
        self._plot_final_PSFs()
        
        # 6. Finalize (Extract OL/DL arrays, compute 1D profiles)
        self._finalize_full_field_results(astIndex)
        
        self.firstSimCall = False

    # ----------------------------------------------------------------------------
    # --- ABSTRACT METHODS (To be implemented by P3 / TipTorch subclasses) ---
    # ----------------------------------------------------------------------------

    @abstractmethod
    def _configure_LO_parameters(self, astIndex):
        """Configure asterism geometry and LO matrices if LO is active."""
        pass

    @abstractmethod
    def _prepare_static_PSF_state(self, astIndex):
        """Initialize the backend model (Fourier/NN) and cache geometry."""
        pass

    @abstractmethod
    def _compute_LO_terms(self, astIndex) -> dict:
        """
        PURE-ISH FUNCTION: Computes Low Order terms based on current geometry.
        Must return a dictionary containing standard keys (e.g., 'HO_res', 'LO_res', 'GF_res', 'GFinPSD').
        """
        pass

    @abstractmethod
    def _generate_final_PSF(self, astIndex, lo_data):
        """
        Generates the final PSFs utilizing the HO state and the provided lo_data.
        Must populate the `self.cubeResultsArray` (standard numpy array).
        """
        pass

    @abstractmethod
    def _finalize_full_field_results(self, astIndex):
        """
        Computes Open-Loop (OL) and Diffraction-Limited (DL) reference PSFs.
        Must populate `self.psf_ol_array` and `self.psf_dl_array` (standard numpy arrays),
        then call `self.computePSF1D()`.
        """
        pass

    @abstractmethod
    def computeMetrics(self):
        """
        Computes Strehl Ratio (SR), FWHM, and Encircled Energy (EE) for the generated PSFs.
        Must populate `self.sr`, `self.fwhm`, and `self.ee`.
        """
        pass

    @abstractmethod
    def _plot_final_PSFs(self):
        """Handles backend-specific rendering (if self.doPlot is True)."""
        pass

    # ----------------------------------------------------------------------------
    # --- CONCRETE METHODS (Universal Logic shared across all backends) ---
    # ----------------------------------------------------------------------------

    def _unpack_LO_data(self, lo_data: dict):
        """Safely unpacks the returned LO dictionary into instance attributes."""
        if not lo_data:
            return
        self.HO_res = lo_data.get('HO_res', self.HO_res)
        self.LO_res = lo_data.get('LO_res', self.LO_res)
        self.GF_res = lo_data.get('GF_res', self.GF_res)
        self.GFinPSD = lo_data.get('GFinPSD', self.GFinPSD)
        self.cov_ellipses = lo_data.get('cov_ellipses', getattr(self, 'cov_ellipses', None))

    def computePSF1D(self):
        """Universal 1D radial profile computation based on standardized numpy arrays."""
        if self.cubeResultsArray is None:
            raise RuntimeError("cubeResultsArray must be populated before calling computePSF1D")
            
        psf1d = []
        psf1d_radius_list_list = []
        use_polar_interp = self.SupSamp and self.SupSamp[1] == 2
        polar_grid, r_vals_interp = None, None
        
        # Precompute polar grid if required
        if use_polar_interp:
            step_interp = self.SupSamp[0]
            # Handle both (nWvl, nPointings, N, N) and (nPointings, N, N)
            first_psf = self.cubeResultsArray[0, 0] if self.nWvl > 1 else self.cubeResultsArray[0]
            center = np.unravel_index(np.argmax(first_psf), first_psf.shape)
            maxradius = self.psInMas * (first_psf.shape[0] / 2)
            r_vals_interp, polar_grid = precompute_polar_grid(step=step_interp,
                                                              pixelscale=self.psInMas,
                                                              maxradius=maxradius,
                                                              center=center)
        for i in range(self.nWvl):
            cubeResults = self.cubeResultsArray[i] if self.nWvl > 1 else self.cubeResultsArray
            psf1dList, psf1d_radius_list = [], []
            
            for psf in cubeResults:
                psfRadius = psf.shape[0] / 2
                center = np.unravel_index(np.argmax(psf), psf.shape)
                rr, radialprofile, _ = radial_profile(psf,
                                                      ext=0,
                                                      pixelscale=self.psInMas,
                                                      ee=True,
                                                      center=center,
                                                      stddev=False,
                                                      binsize=None,
                                                      maxradius=self.psInMas * psfRadius,
                                                      normalize='total',
                                                      nargout=2,
                                                      supersamp=self.SupSamp, 
                                                      polar_grid=polar_grid,
                                                      r_vals=r_vals_interp,
                                                      verbose=self.verbose)
                psf1dList.append(radialprofile)
                psf1d_radius_list.append(rr)
                
            psf1d.append(psf1dList)
            psf1d_radius_list_list.append(psf1d_radius_list)
            
        self.psf1d = np.asarray(psf1d)
        self.psf1d_radius = np.asarray(psf1d_radius_list_list[0][0])
        # All PSFs share the same shape and psInMas, so radii must be identical.
        assert all(
            np.allclose(psf1d_radius_list_list[i][j], self.psf1d_radius)
            for i in range(self.nWvl)
            for j in range(len(psf1d_radius_list_list[i]))
        ), "Radial profile radii differ across PSFs — check PSF shapes and pixel scale consistency."
        self.psf1d_data = np.vstack((np.asarray(psf1d_radius_list_list), self.psf1d))

    def savePSFprofileJSON(self):
        """Universal JSON serialization for PSF profiles."""
        now = datetime.now()
        psf_data = {
            'radius': self.psf1d_radius.tolist(),
            'psf': self.psf1d.tolist()
        }
        filename = os.path.join(self.outputDir, self.outputFile + '1D_PSF.json')
        jsondict = {
            'execution_infos': {
                'TIME': now.strftime("%Y%m%d_%H%M%S"),
                'TIPTOP version': __version__
            },
            'infos': self.my_data_map,
            'psf': psf_data
        }
        with open(filename, 'w') as f:
            json.dump(jsondict, f)

    def saveResults(self):
        """
        Universal FITS saving method. 
        It strictly requires the backend to have populated the standardized NumPy arrays
        (self.cubeResultsArray, self.psf_ol_array, self.psf_dl_array, self.psf1d_data)
        ensuring complete decoupling from specific tensors or Field objects.
        """
        if any(v is None for v in [self.cubeResultsArray, self.psf_ol_array, self.psf_dl_array, self.psf1d_data]):
            raise RuntimeError("Backend failed to populate standard array containers before saving.")

        hdul1 = fits.HDUList()
        hdul1.append(fits.PrimaryHDU())
        hdul1.append(fits.ImageHDU(data=self.cubeResultsArray))
        hdul1.append(fits.ImageHDU(data=self.psf_ol_array))
        hdul1.append(fits.ImageHDU(data=self.psf_dl_array))
        if self.savePSDs and self.psd_array is not None:
            hdul1.append(fits.ImageHDU(data=self.psd_array))
        hdul1.append(fits.ImageHDU(data=self.psf1d_data))

        now = datetime.now()
        hdr0 = hdul1[0].header       
        hdr0['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        hdr0['TIPTOP_V'] = __version__
        
        # Extract configuration metadata
        for key_primary, section in self.my_data_map.items():
            for key_secondary, temp in section.items():
                if isinstance(temp, list):
                    for iii, elem in enumerate(temp):
                        if isinstance(elem, list):
                            for jjj, elem2 in enumerate(elem):
                                add_hdr_keyword(hdr0, key_primary, key_secondary, elem2, iii=str(iii), jjj=str(jjj))
                        else:                        
                            add_hdr_keyword(hdr0, key_primary, key_secondary, elem, iii=str(iii))
                else:
                    add_hdr_keyword(hdr0, key_primary, key_secondary, temp)

        # Main Data Header
        hdr1 = hdul1[1].header
        hdr1['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        hdr1['CONTENT'] = "PSF CUBE"
        hdr1['SIZE'] = str(self.cubeResultsArray.shape)
        
        if self.nWvl > 1:
            for i in range(self.nWvl):
                hdr1[f'WL_NM{str(i).zfill(3)}'] = str(int(self.wvl[i] * 1e9))
        else:
            hdr1['WL_NM'] = str(int(self.wvl[0] * 1e9))
            
        hdr1['PIX_MAS'] = str(self.psInMas)
        hdr1['CC'] = f"CARTESIAN COORD. IN ASEC OF THE {self.pointings.shape[1]} SOURCES"
        for i in range(self.pointings.shape[1]):
            hdr1[f'CCX{str(i).zfill(4)}'] = float(np.round(self.pointings[0, i], 3).item())
            hdr1[f'CCY{str(i).zfill(4)}'] = float(np.round(self.pointings[1, i], 3).item())
            
        # SANITIZATION: Strict float casting to avoid FITS Header TypeErrors
        if self.HO_res is not None:
            hdr1['RESH'] = "High Order residual in nm RMS"
            # Handling vectors of residuals
            ho_array = np.atleast_1d(self.HO_res)
            for i in range(ho_array.shape[0]):
                hdr1[f'RESH{str(i).zfill(4)}'] = float(np.round(ho_array[i].item(), 3))
                
        if self.LO_res is not None:
            hdr1['RESL'] = "Low Order residual in nm RMS"
            lo_array = np.atleast_1d(self.LO_res)
            for i in range(lo_array.shape[0]):
                hdr1[f'RESL{str(i).zfill(4)}'] = float(np.round(lo_array[i].item(), 3))
                
        if self.GF_res is not None:
            hdr1['RESF'] = "Global Focus residual in nm RMS (included in PSD)"
            hdr1['RESF0000'] = float(np.round(np.asarray(self.GF_res).item(), 3))
            
        # Add Metrics to Header
        if self.addSrAndFwhm:
            # Ensure computeMetrics has populated sr, fwhm, ee
            if not self.sr or not self.fwhm:
                self.computeMetrics()
                
            for i in range(self.nWvl):
                wTxt, fTxt, eTxt, Nfill = (f'W{str(i).zfill(2)}', 'FW', 'EE', 2) if self.nWvl > 1 \
                                          else ('', 'FWHM', f'EE{int(self.eeRadiusInMas)}', 4)
                
                cube_slice = self.cubeResultsArray[i] if self.nWvl > 1 else self.cubeResultsArray
                sr_slice = self.sr[i] if self.nWvl > 1 else self.sr
                fwhm_slice = self.fwhm[i] if self.nWvl > 1 else self.fwhm
                ee_slice = self.ee[i] if self.nWvl > 1 else self.ee

                for j in range(cube_slice.shape[0]):
                    sr_val = float(np.round(np.asarray(cpuArray(sr_slice[j])).item(), 5))
                    fwhm_val = float(np.round(np.asarray(cpuArray(fwhm_slice[j])).item(), 3))
                    ee_val = float(np.round(np.asarray(cpuArray(ee_slice[j])).item(), 5))

                    hdr1[f'SR{str(j).zfill(Nfill)}{wTxt}'] = sr_val
                    hdr1[f'{fTxt}{str(j).zfill(Nfill)}{wTxt}'] = fwhm_val
                    hdr1[f'{eTxt}{str(j).zfill(Nfill)}{wTxt}'] = ee_val

        # Secondary Headers
        hdul1[2].header.update({'TIME': now.strftime("%Y%m%d_%H%M%S"),
                                'CONTENT': "OPEN-LOOP PSF",
                                'SIZE': str(self.psf_ol_array.shape)})
        hdul1[3].header.update({'TIME': now.strftime("%Y%m%d_%H%M%S"),
                                'CONTENT': "DIFFRACTION LIMITED PSF",
                                'SIZE': str(self.psf_dl_array.shape)})
        
        idx_offset = 4
        if self.savePSDs and self.psd_array is not None:
            hdul1[4].header.update({'TIME': now.strftime("%Y%m%d_%H%M%S"),
                                    'CONTENT': "High Order PSD",
                                    'SIZE': str(self.psd_array.shape)})
            idx_offset = 5
            
        hdul1[idx_offset].header.update({'TIME': now.strftime("%Y%m%d_%H%M%S"),
                                        'CONTENT': "Final PSFs profiles",
                                        'SIZE': str(self.psf1d_data.shape)})
        if self.SupSamp:
            hdul1[idx_offset].header['SAMP_MAS'] = str(self.SupSamp[0])

        hdul1.writeto(os.path.join(self.outputDir, self.outputFile + '.fits'), overwrite=True)
        
        if self.verbose:
            print("Output cube shape:", self.cubeResultsArray.shape)
            print("Output dtype:", self.cubeResultsArray.dtype)
