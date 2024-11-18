from p3.aoSystem.fourierModel import *
from p3.aoSystem.FourierUtils import *
from mastsel import *

from .tiptopUtils import *

from matplotlib import cm
import matplotlib as mpl
norm = mpl.colors.Normalize(vmin=0, vmax=1)
rc("text", usetex=False)

class baseSimulation(object):
    
    def raiseMissingRequiredOpt(self,sec,opt):
        raise ValueError("'{}' is missing from section '{}'"
                         .format(opt,sec))
        
    def raiseMissingRequiredSec(self,sec):
        raise ValueError("The section '{}' is missing from the parameter file"
                         .format(sec))
        
    def raiseNotSameLength(self,sec,opt):
        raise ValueError("'{}' in section '{}' must have the same length"
                         .format(opt,sec))
    
    def check_section_key(self, primary):        
        return primary in self.my_data_map.keys()
    
    def check_config_key(self, primary, secondary):
        if primary in self.my_data_map.keys():
            return secondary in self.my_data_map[primary].keys()
        else:
            return False
    
    def __init__(self, path, parametersFile, outputDir, outputFile, doConvolve=True,
                          doPlot=False, addSrAndFwhm=False,
                          verbose=False, getHoErrorBreakDown=False,
                          savePSDs=False, ensquaredEnergy=False,
                          eeRadiusInMas=50):
        self.firstSimCall =True
        if verbose: np.set_printoptions(precision=3)
        self.doConvolveAsterism = True
        self.pointings_FWHM_mas = None
        # copy the parameters in state vars
        self.path = path
        self.parametersFile = parametersFile
        self.outputDir = outputDir
        self.outputFile = outputFile
        self.doPlot = doPlot
        self.doConvolve = doConvolve
        self.addSrAndFwhm = addSrAndFwhm
        self.verbose = verbose
        self.getHoErrorBreakDown = getHoErrorBreakDown
        self.savePSDs = savePSDs
        self.ensquaredEnergy = ensquaredEnergy
        self.eeRadiusInMas = eeRadiusInMas
#        if self.returnRes and self.doPlot:
#            print('WARNING: returnRes and doPlot cannot both be True, setting doPlot to False.')
#            self.doPlot = False

        # get the system description (stored in my_data_map) from the ini/yml file
        fullPathFilename_ini = os.path.join(self.path, self.parametersFile + '.ini')
        fullPathFilename_yml = os.path.join(self.path, self.parametersFile + '.yml')
        if os.path.exists(fullPathFilename_yml):
            self.fullPathFilename = fullPathFilename_yml
            with open(fullPathFilename_yml) as f:
                my_yaml_dict = yaml.safe_load(f)        
            self.my_data_map = my_yaml_dict
        elif os.path.exists(fullPathFilename_ini):
            self.fullPathFilename = fullPathFilename_ini
            config = ConfigParser()
            config.optionxform = str
            config.read(fullPathFilename_ini)
            self.my_data_map = {} 
            for section in config.sections():
                self.my_data_map[section] = {}
                for name,value in config.items(section):
                    self.my_data_map[section].update({name:eval(value)})
            
            #Verify the presence of parameters called in TIPTOP before they are verified 
            #in P3 or in MASTSEL
            if not self.check_section_key('telescope'):
                self.raiseMissingRequiredSec('telescope')
                
            if not self.check_config_key('telescope','TelescopeDiameter'):
                self.raiseMissingRequiredOpt('telescope','TelescopeDiameter')
            
            if not self.check_config_key('telescope','glFocusOnNGS'):
                self.my_data_map['telescope']['glFocusOnNGS'] = False
            
            if not self.check_section_key('sources_science') :
                self.raiseMissingRequiredSec('sources_science') 
            
            if not self.check_config_key('sources_science','Wavelength'):
                self.raiseMissingRequiredOpt('sources_science', 'Wavelength')
            
            if not self.check_config_key('sources_science','Zenith'):
                #In P3.aoSystem this is optionnal, to remain consistent it is optionnal here too
                self.my_data_map['sources_science']['Zenith'] = [0.0]
            
            if not self.check_config_key('sources_science','Azimuth'):
                #In P3.aoSystem this is optionnal, to remain consistent it is optionnal here too
                self.my_data_map['sources_science']['Azimuth'] = [0.0]
            
            if (len(self.my_data_map['sources_science']['Zenith']) != 
                len(self.my_data_map['sources_science']['Azimuth'])):
                self.raiseNotSameLength('sources_science', ['Zenith','Azimuth'])
            
            #TODO should an error be raised if sensor_LO is defined but not source_LO or vice versa?
            if self.check_section_key('sources_LO') and not self.check_section_key('sensor_LO'):
                raise KeyError("'sensor_LO' must be defined if 'sources_LO' is defined.")
            elif not self.check_section_key('sources_LO') and self.check_section_key('sensor_LO'):
                raise KeyError("'sources_LO' must be defined if 'sensor_LO' is defined.")
            #If both are defined we can proceed.
            elif self.check_section_key('sources_LO') and self.check_section_key('sensor_LO'):
                if not self.check_config_key('sources_LO', 'Wavelength'):
                    self.raiseMissingRequiredOpt('sources_LO', 'Wavelength')
                
                if not self.check_config_key('sources_LO','Zenith'):
                    self.my_data_map['sources_LO']['Zenith'] = [0.0]
                    
                if not self.check_config_key('sources_LO','Azimuth'):
                    self.my_data_map['sources_LO']['Azimuth'] = [0.0]
                
                if (len(self.my_data_map['sources_LO']['Zenith']) != 
                    len(self.my_data_map['sources_LO']['Azimuth'])):
                    self.raiseNotSameLength('sources_LO', ['Zenith','Azimuth'])
                
                if not self.check_config_key('sensor_LO', 'NumberPhotons'):
                    raise self.raiseMissingRequiredOpt('sensor_LO', 'NumberPhotons')
                
                if not self.check_section_key('RTC'):
                    self.raiseMissingRequiredSec('RTC')
                elif not self.check_config_key('RTC', 'SensorFrameRate_LO'):
                    self.raiseMissingRequiredOpt('RTC', 'SensorFrameRate_LO')
            
        else:
            raise FileNotFoundError('No .yml or .ini can be found in '+ self.path)

        self.tel_radius = self.my_data_map['telescope']['TelescopeDiameter']/2  # mas
        wvl_temp = self.my_data_map['sources_science']['Wavelength']
        if isinstance(wvl_temp, list):
            self.wvl = wvl_temp[0]  # lambda
            if len(wvl_temp) > 1:
                print('WARNING: TIPTOP is not supporting multi-wavelength case, it will provide \
                       PSF for the first wavelength of the list',wvl_temp,'.')
        else:
            self.wvl = wvl_temp     # lambda
        self.zenithSrc  = self.my_data_map['sources_science']['Zenith']
        self.azimuthSrc = self.my_data_map['sources_science']['Azimuth']
        self.pointings = polarToCartesian(np.array( [self.zenithSrc, self.azimuthSrc]))
        self.xxSciencePointigs         = self.pointings[0,:]
        self.yySciencePointigs         = self.pointings[1,:]
        # it checks if LO parameters are set and then it acts accordingly
        if 'sensor_LO' in self.my_data_map.keys():
            self.LOisOn = True
            if self.verbose: print('LO part is present')
        else:
            self.LOisOn = False
            self.nNaturalGS_field = 0
            if self.verbose: print('LO part is not present')
            
        # initialize self.jitter_FWHM variable with a default value
        self.jitter_FWHM = None
        if 'jitter_FWHM' in self.my_data_map['telescope'].keys():
            self.jitter_FWHM = self.my_data_map['telescope']['jitter_FWHM']  
            
        self.addFocusError = self.my_data_map['telescope']['glFocusOnNGS']
        self.GFinPSD = False
        if (not self.check_section_key('sensor_Focus')) and self.addFocusError and max(self.my_data_map['sensor_LO']['NumberLenslets']) == 1:
            raise ValueError("[telescope] glFocusOnNGS (that is focus correction with NGS) is available only if NGS/Focus WFSs have more than one sub-aperture")


    def configLO(self, astIndex=None):
        self.cartSciencePointingCoords = np.dstack( (self.xxSciencePointigs, self.yySciencePointigs) ).reshape(-1, 2)
        # Here we assume the same wavelenght for all the phon counts of the stars in the asterism
        LO_wvl_temp = self.my_data_map['sources_LO']['Wavelength']

        if isinstance(LO_wvl_temp, list):
            self.LO_wvl = LO_wvl_temp[0]  # lambda
        else:
            self.LO_wvl = LO_wvl_temp     # lambda

        self.LO_psInMas      = self.my_data_map['sensor_LO']['PixelScale']
        self.LO_zen_field    = self.my_data_map['sources_LO']['Zenith']
        self.LO_az_field     = self.my_data_map['sources_LO']['Azimuth']
        self.LO_fluxes_field = self.my_data_map['sensor_LO']['NumberPhotons']
        self.LO_freqs_field  = self.my_data_map['RTC']['SensorFrameRate_LO']
        if not isinstance(self.LO_freqs_field, list):
            self.LO_freqs_field  = [self.LO_freqs_field] * len(self.LO_zen_field)

        if self.check_section_key('sensor_Focus'):
            self.Focus_fluxes4s_field = self.my_data_map['sensor_Focus']['NumberPhotons']
            self.Focus_psInMas         = self.my_data_map['sensor_Focus']['PixelScale']
            Focus_wvl_temp = self.my_data_map['sources_LO']['Wavelength']
            if isinstance(Focus_wvl_temp, list):
                self.Focus_wvl = Focus_wvl_temp[0]  # lambda
            else:
                self.Focus_wvl = Focus_wvl_temp     # lambda
        else:
            self.Focus_fluxes4s_field = self.LO_fluxes_field
            self.Focus_psInMas         = self.LO_psInMas
            self.Focus_wvl             = self.LO_wvl
        if self.check_config_key('RTC','SensorFrameRate_Focus'):
            self.Focus_freqs_field  = self.my_data_map['RTC']['SensorFrameRate_Focus']
            if not isinstance(self.Focus_freqs_field, list):
                self.Focus_freqs_field  = [self.Focus_freqs_field] * len(self.LO_zen_field)
        else:
            self.Focus_freqs_field = self.LO_freqs_field

        self.NGS_fluxes_field = []
        polarNGSCoordsList = []
        for aFr, aFlux, aZen, aAz in zip(self.LO_freqs_field, self.LO_fluxes_field, self.LO_zen_field, self.LO_az_field):
            polarNGSCoordsList.append([aZen, aAz])
            self.NGS_fluxes_field.append(aFlux*aFr)
        self.Focus_fluxes_field = []
        for aFrF, aFluxF in zip(self.Focus_freqs_field, self.Focus_fluxes4s_field):
            self.Focus_fluxes_field.append(aFluxF*aFrF)
        polarNGSCoords     = np.asarray(polarNGSCoordsList)
        self.nNaturalGS_field  = len(self.LO_zen_field)
        cartNGSCoordsList = []
        for i in range(self.nNaturalGS_field):
            cartNGSCoordsList.append(polarToCartesian(polarNGSCoords[i,:]))

        self.cartNGSCoords_field = np.asarray(cartNGSCoordsList)            
        self.currentAsterismIndices = list(range(len(self.LO_zen_field)))
        self.setAsterismData()


    def setAsterismData(self):
        self.LO_zen_asterism = []
        for iid in self.currentAsterismIndices:
            self.LO_zen_asterism.append(self.LO_zen_field[iid])
        self.LO_az_asterism = []
        for iid in self.currentAsterismIndices:
            self.LO_az_asterism.append(self.LO_az_field[iid])
        self.LO_fluxes_asterism = []
        for iid in self.currentAsterismIndices:
            self.LO_fluxes_asterism.append(self.LO_fluxes_field[iid])
        self.LO_freqs_asterism = []
        for iid in self.currentAsterismIndices:
            self.LO_freqs_asterism.append(self.LO_freqs_field[iid])
        self.NGS_fluxes_asterism = []
        for iid in self.currentAsterismIndices:
            self.NGS_fluxes_asterism.append(self.NGS_fluxes_field[iid])
        self.Focus_fluxes_asterism = []
        for iid in self.currentAsterismIndices:
            self.Focus_fluxes_asterism.append(self.Focus_fluxes_field[iid])
        self.cartNGSCoords_asterism = []
        for iid in self.currentAsterismIndices:
            self.cartNGSCoords_asterism.append(self.cartNGSCoords_field[iid])


    def saveResults(self):
        # save PSF cube in fits
        hdul1 = fits.HDUList()
        hdul1.append(fits.PrimaryHDU())        
        hdul1.append(fits.ImageHDU(data=self.cubeResultsArray))
        hdul1.append(fits.ImageHDU(data=self.psfOL.sampling)) # append open-loop PSF
        hdul1.append(fits.ImageHDU(data=self.psfDL.sampling)) # append diffraction limited PSF
        if self.savePSDs:
            hdul1.append(fits.ImageHDU(data=cpuArray(self.PSD))) # append high order PSD
        # header
        hdr0 = hdul1[0].header
        now = datetime.now()
        hdr0['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        # parameters in the header
        for key_primary in self.my_data_map:
            for key_secondary in self.my_data_map[key_primary]:
                temp = self.my_data_map[key_primary][key_secondary]
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
        hdr1['WL_NM'] = str(int(self.wvl*1e9))
        hdr1['PIX_MAS'] = str(self.psInMas[0])
        hdr1['CC'] = "CARTESIAN COORD. IN ASEC OF THE "+str(self.pointings.shape[1])+" SOURCES"
        for i in range(self.pointings.shape[1]):
            hdr1['CCX'+str(i).zfill(4)] = self.pointings[0,i]
            hdr1['CCY'+str(i).zfill(4)] = self.pointings[1,i]
        if hasattr(self,'HO_res'):
            hdr1['RESH'] = "High Order residual in nm RMS"
            for i in range(self.HO_res.shape[0]):
                hdr1['RESH'+str(i).zfill(4)] = str(self.HO_res[i])
        if hasattr(self,'LO_res'):
            hdr1['RESL'] = "Low Order residual in nm RMS"
            for i in range(self.LO_res.shape[0]):
                hdr1['RESL'+str(i).zfill(4)] = str(self.LO_res[i])
        if hasattr(self,'GF_res'):
            hdr1['RESF'] = "Global Focus residual in nm RMS (included in PSD)"
            hdr1['RESF0000'] = str(self.GF_res)
        if self.addSrAndFwhm:
            for i in range(self.cubeResultsArray.shape[0]):
                hdr1['SR'+str(i).zfill(4)]   = float(getStrehl(self.cubeResultsArray[i,:,:], self.fao.ao.tel.pupil, self.fao.freq.sampRef, method='otf'))
            for i in range(self.cubeResultsArray.shape[0]):
                hdr1['FWHM'+str(i).zfill(4)] = getFWHM(self.cubeResultsArray[i,:,:], self.psInMas[0], method='contour', nargout=1)
            for i in range(self.cubeResultsArray.shape[0]):
                if self.ensquaredEnergy:
                    ee = cpuArray(getEnsquaredEnergy(self.cubeResultsArray[i,:,:]))
                    rr = np.arange(1, ee.shape[0]*2, 2) * self.psInMas[0] * 0.5
                else:
                    ee,rr = getEncircledEnergy(self.cubeResultsArray[i,:,:], pixelscale=self.psInMas[0], center=(self.fao.ao.cam.fovInPix/2,self.fao.ao.cam.fovInPix/2), nargout=2)
                ee_at_radius_fn = interp1d(rr, ee, kind='cubic', bounds_error=False)
                hdr1['EE'+str(np.round(self.eeRadiusInMas))+str(i).zfill(4)] = ee_at_radius_fn(self.eeRadiusInMas).take(0)

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

        if self.savePSDs:
            # header of the PSD
            hdr4 = hdul1[4].header
            hdr4['TIME'] = now.strftime("%Y%m%d_%H%M%S")
            hdr4['CONTENT'] = "High Order PSD"
            hdr4['SIZE'] = str(self.PSD.shape)

        hdul1.writeto( os.path.join(self.outputDir, self.outputFile + '.fits'), overwrite=True)
        if self.verbose:
            print("Output cube shape:", self.cubeResultsArray.shape)


    def psdSetToPsfSet(self, N, freq_range, dk, mask, inputPSDs, wavelength, pixelscale, npixel, scaleFactor=1, oversampling=1):
        sources_SR = []
        psdArray = []
        psfLongExpArr = []
        sources_FWHM_mas = []

        # Get the Telecope plus Static WFE OTF if defined in self
        if hasattr(self, "otfStatic"):
            otf_tel = self.otfStatic[0]
            # above only 0 is used because this method does not yet support multi wavelength PSF generation
        else:
            otf_tel = None

        i = 0
        for computedPSD in inputPSDs:
            # Get the PSD at the NGSs positions at the sensing wavelength
            # computed PSD from fao are given in nm^2, i.e they are multiplied by dk**2 already
            psd          = Field(wavelength, N, freq_range, 'rad')
            psd.sampling = computedPSD / dk**2 # the PSD must be provided in m^2.m^2
            psdArray.append(psd)
            # Get the PSF
            if isinstance(mask, list):
                psfLE = longExposurePsf(mask[i], psd, otf_tel = otf_tel )
            else:
                psfLE = longExposurePsf(mask, psd, otf_tel = otf_tel )
            
            # It rebins the PSF if oversampling is greater than 1
            if oversampling > 1:
                temp = np.array(psfLE.sampling)
                nTemp = int(oversampling)
                nOut = int(temp.shape[0]/nTemp)
                psfLE.sampling = temp.reshape((nOut,nTemp,nOut,nTemp)).mean(3).mean(1)
            # It cuts the PSF if the PSF is larger than the requested dimension
            if psfLE.sampling.shape[0] > npixel:
                psfLE.sampling = psfLE.sampling[int(psfLE.sampling.shape[0]/2-npixel/2):int(psfLE.sampling.shape[0]/2+npixel/2),
                                                int(psfLE.sampling.shape[1]/2-npixel/2):int(psfLE.sampling.shape[1]/2+npixel/2)]
            psfLongExpArr.append(psfLE)
            # Get SR and FWHM in mas at the NGSs positions at the sensing wavelength
            s1           = cpuArray(computedPSD).sum()
            SR           = np.exp(-s1*scaleFactor) # Strehl-ratio at the sensing wavelength

            sources_SR.append(SR)
            FWHMx,FWHMy  = getFWHM( psfLE.sampling, pixelscale, method='contour', nargout=2)
            FWHM         = np.sqrt(FWHMx*FWHMy) #max(FWHMx, FWHMy) #0.5*(FWHMx+FWHMy) #average over major and minor axes

            # note : the uncertainities on the FWHM seems to create a bug in mavisLO
            sources_FWHM_mas.append(FWHM)
            if self.verbose:
                print('SR(@',int(wavelength*1e9),'nm)        :', "%.5f" % SR)
                print('FWHM(@',int(wavelength*1e9),'nm) [mas]:', "%.3f" % FWHM)
            i += 1

        return sources_SR, psdArray, psfLongExpArr, sources_FWHM_mas


    def computeOL_PSD(self):
        # OPEN-LOOP PSD
        k   = np.sqrt(self.fao.freq.k2_)
        pf  = FourierUtils.pistonFilter(self.fao.ao.tel.D,k)
        psdOL = Field(self.wvl, self.N, self.freq_range, 'rad')
        temp = self.fao.ao.atm.spectrum(k) * pf
        psdOL.sampling = arrayP3toMastsel(self.fao.ao.atm.spectrum(k) * pf * (self.fao.freq.wvlRef/np.pi)**2) # the PSD must be provided in m^2.m^2
        # Get the OPEN-LOOP PSF
        self.psfOL = longExposurePsf(self.mask, psdOL)
        # It cuts the PSF if the PSF is larger than the requested dimension (N>nPixPSF)
#        if self.psfOL.sampling.shape[0] > nPixPSF:
#            self.psfOL.sampling = self.psfOL.sampling[int(self.psfOL.sampling.shape[0]/2-nPixPSF/2):int(self.psfOL.sampling.shape[0]/2+nPixPSF/2),
#                                            int(self.psfOL.sampling.shape[1]/2-nPixPSF/2):int(self.psfOL.sampling.shape[1]/2+nPixPSF/2)]
#        self.psfOL.sampling = cpuArray(psdOL.sampling)
        if self.doPlot:
            fig, ax1 = plt.subplots(1,1)
            im = ax1.imshow(np.log(np.abs(self.psfOL.sampling) + 1e-20), cmap='hot')
            ax1.set_title('open loop PSF', color='black')


    def computeDL_PSD(self):
        # DIFFRACTION LIMITED PSD
        psdDL = Field(self.wvl, self.N, self.freq_range, 'rad')
        self.psfDL = longExposurePsf(self.mask, psdDL)
        # It cuts the PSF if the PSF is larger than the requested dimension (N>nPixPSF)
        if self.psfDL.sampling.shape[0] > self.nPixPSF:
            self.psfDL.sampling = self.psfDL.sampling[int(self.psfDL.sampling.shape[0]/2-self.nPixPSF/2):int(self.psfDL.sampling.shape[0]/2+self.nPixPSF/2),
                                            int(self.psfDL.sampling.shape[1]/2-self.nPixPSF/2):int(self.psfDL.sampling.shape[1]/2+self.nPixPSF/2)]
        if self.doPlot:
            fig, ax2 = plt.subplots(1,1)
            im = ax2.imshow(np.log(np.abs(self.psfDL.sampling) + 1e-20), cmap='hot')
            ax2.set_title('diffraction limited PSF', color='black')


    def finalConvolution(self):
        self.cov_ellipses = self.mLO.ellipsesFromCovMats(self.Ctot)
        if self.verbose:
            for n in range(self.cov_ellipses.shape[0]):
                print('cov_ellipses #',n,': ', self.cov_ellipses[n,:], ' (unit: rad, mas, mas)')
            print('******** FINAL CONVOLUTION')
        # FINAl CONVOLUTION
        # Optimization: Non NEED to perform this convolutions if this is the asterism selection procedure ???
        for ellp, psfLongExp in zip(self.cov_ellipses, self.psfLongExpPointingsArr):
            resSpec = residualToSpectrum(ellp, self.wvl, self.nPixPSF, 1/(self.fao.ao.cam.fovInPix * self.psInMas[0]))
            temp = convolve(psfLongExp, resSpec)
            if self.jitter_FWHM is not None:
                if isinstance(self.jitter_FWHM, list):
                    ellpJ = [self.jitter_FWHM[2], sigma_from_FWHM(self.jitter_FWHM[0]), sigma_from_FWHM(self.jitter_FWHM[1])]
                else:
                    ellpJ = [0, sigma_from_FWHM(self.jitter_FWHM), sigma_from_FWHM(self.jitter_FWHM)]
                resSpecJ = residualToSpectrum(ellpJ, self.wvl, self.nPixPSF, 1/(self.fao.ao.cam.fovInPix * self.psInMas[0]))
                temp = convolve(temp, resSpecJ)
            self.results.append(temp)


    def finalPSF(self,astIndex):
        if astIndex is None or self.firstSimCall:
            # ----------------------------------------------------------------------------
            ## HO PSF
            if self.verbose:
                print('******** HO PSF')
            pointings_SR, psdPointingsArray, psfLongExpPointingsArr, self.pointings_FWHM_mas = self.psdSetToPsfSet(self.N,
                                                                                                         self.freq_range,
                                                                                                         self.dk,
                                                                                                         self.mask,
                                                                                                         arrayP3toMastsel(self.PSD[0:self.nPointings]),
                                                                                                         self.wvl,
                                                                                                         self.psInMas[0],
                                                                                                         self.nPixPSF,
                                                                                                         scaleFactor=(2*np.pi*1e-9/self.wvl)**2,
                                                                                                         oversampling=self.oversampling)
            self.psfLongExpPointingsArr = psfLongExpPointingsArr

            # ----------------------------------------------------------------------------
            ## computation of the HO error (this is fixed for the simulation)
            self.HO_res = np.sqrt(np.sum(self.PSD[0:self.nPointings],axis=(1,2)))

        # ------------------------------------------------------------------------
        ## final PSFs computation after optional convolution with jitter kernels
        if self.LOisOn:
            if self.doConvolve:
                if self.doConvolveAsterism:
                    self.finalConvolution()
                else:
                    self.cov_ellipses = self.mLO.ellipsesFromCovMats(self.Ctot)
            else:
                for psfLongExp in self.psfLongExpPointingsArr:
                    self.results.append(psfLongExp)
        else:
            for psfLongExp in self.psfLongExpPointingsArr:
                if self.jitter_FWHM is not None:
                    if isinstance(self.jitter_FWHM, list):
                        ellp = [self.jitter_FWHM[2], sigma_from_FWHM(self.jitter_FWHM[0]), sigma_from_FWHM(self.jitter_FWHM[1])]
                    else:
                        ellp = [0, sigma_from_FWHM(self.jitter_FWHM), sigma_from_FWHM(self.jitter_FWHM)]
                    self.results.append(convolve(psfLongExp,
                                   residualToSpectrum(ellp, self.wvl, self.nPixPSF, 1/(self.fao.ao.cam.fovInPix * self.psInMas[0]))))
                else:
                    self.results.append(psfLongExp)


    def ngsPSF(self):
        # -----------------------------------------------------------------
        # PSD and sub-aperture mask for NGS directions
        PSD_NGS = arrayP3toMastsel(self.PSD[-self.nNaturalGS_field:])
        k  = np.sqrt(self.fao.freq.k2_)
        # Define the LO sub-aperture shape
        nSA = self.my_data_map['sensor_LO']['NumberLenslets']
        pupilSidePix = int(self.fao.ao.tel.pupil.shape[0])
        saMask = np.zeros((pupilSidePix,pupilSidePix))
        if len(nSA) == self.nNaturalGS_field:
            self.maskLO = []
            nMaskLO = self.nNaturalGS_field
        else:
            nMaskLO = 1
        for i in range(nMaskLO):
            saSideM = 2*self.tel_radius/nSA[i]
            saSidePix = int(pupilSidePix/nSA[i])
            if nSA[i] == 1:
                # LO mask
                if nMaskLO > 1:
                    self.maskLO.append(self.mask)
                else:
                    self.maskLO = self.mask
            else:
                # piston filter for the sub-aperture size
                pf = FourierUtils.pistonFilter(self.fao.ao.tel.D/nSA[i],k)
                PSD_NGS[i] = PSD_NGS[i] * pf
                # LO mask
                maskLO = Field(self.wvl, self.N, self.grid_diameter)
                if nSA[i] == 2:
                    saMask[0:saSidePix,0:saSidePix] = 1
                    saMask *= cpuArray(self.fao.ao.tel.pupil)
                elif nSA[i] == 3:
                    saMask[0:saSidePix,int(pupilSidePix/2-saSidePix/2):int(pupilSidePix/2+saSidePix/2)] = 1
                    saMask *= cpuArray(self.fao.ao.tel.pupil)
                else:
                    saMask[int(pupilSidePix/2-saSidePix/2):int(pupilSidePix/2+saSidePix/2),\
                           int(pupilSidePix/2-saSidePix/2):int(pupilSidePix/2+saSidePix/2)] = 1
                if gpuMastsel:
                    maskLO.sampling = congrid(cp.asarray(saMask), [self.sx, self.sx])
                else:
                    maskLO.sampling = congrid(saMask, [self.sx, self.sx])
                maskLO.sampling = zeroPad(maskLO.sampling, (self.N-self.sx)//2)
                if nMaskLO > 1:
                    self.maskLO.append(maskLO)
                else:
                    self.maskLO = maskLO

        # -----------------------------------------------------------------
        # PSF for NGS directions

        if self.verbose:
            print('******** LO PSF - NGS directions (1 sub-aperture)')

        NGS_SR, psdArray, psfLE_NGS, NGS_FWHM_mas = self.psdSetToPsfSet(self.N,
                                                                   self.freq_range,
                                                                   self.dk,
                                                                   self.maskLO,
                                                                   PSD_NGS,
                                                                   self.LO_wvl,
                                                                   self.psInMas[0],
                                                                   self.nPixPSF,
                                                                   scaleFactor=(2*np.pi*1e-9/self.LO_wvl)**2,
                                                                   oversampling=self.oversampling)

        # -----------------------------------------------------------------
        # Merit functions
        self.NGS_SR_field           = NGS_SR
        self.NGS_FWHM_mas_field     = NGS_FWHM_mas
        self.NGS_EE_field           = []
        idx = 0
        for img in psfLE_NGS:
            ee_,rr_ = getEncircledEnergy(img.sampling, pixelscale=self.psInMas[0],
                                         center=(self.fao.ao.cam.fovInPix/2,self.fao.ao.cam.fovInPix/2), nargout=2)
            ee_ *= 1/np.max(ee_)
            ee_at_radius_fn = interp1d(rr_, ee_, kind='cubic', bounds_error=False)
            # max is used to compute EE on at least a radius of one pixel
            ee_NGS = ee_at_radius_fn(max([NGS_FWHM_mas[idx],self.LO_psInMas]))
            self.NGS_EE_field.append(ee_NGS)
            if self.verbose:
                print('EE                   :', "%.5f" % ee_NGS)
            idx += 1

        # -----------------------------------------------------------------
        # optional Focus error
        if self.addFocusError:
            if 'sensor_Focus' in self.my_data_map.keys():
                if self.verbose:
                    print('Focus sensor is set: computing new PSFs.')
                # -----------------------------------------------------------------
                ## PSD and sub-aperture mask for NGS directions
                nSAfocus = self.my_data_map['sensor_Focus']['NumberLenslets']
                PSD_Focus = arrayP3toMastsel(self.PSD[-self.nNaturalGS_field:])
                if len(nSAfocus) == self.nNaturalGS_field:
                    self.maskFocus = []
                    nMaskFocus = self.nNaturalGS_field
                else:
                    nMaskFocus = 1
                for i in range(nMaskFocus):
                    saSideM = 2*self.tel_radius/nSAfocus[i]
                    saSidePix = int(pupilSidePix/nSAfocus[i])
                    if nSAfocus[i] == 1:
                        # LO mask
                        if nMaskFocus > 1:
                            self.maskFocus.append(self.mask)
                        else:
                            self.maskFocus = self.mask
                    else:
                        ## -----------------------------------------------------------------
                        # --- piston filter for the sub-aperture size
                        pf = FourierUtils.pistonFilter(self.fao.ao.tel.D/nSAfocus[i],k)
                        PSD_Focus[i] = PSD_Focus[i] * pf
                        ## -----------------------------------------------------------------
                        # Focus mask
                        maskFocus = Field(self.wvl, self.N, self.grid_diameter)
                        if nSAfocus[i] == 2:
                            saMask[0:saSidePix,0:saSidePix] = 1
                            saMask *= cpuArray(self.fao.ao.tel.pupil)
                        elif nSAfocus[i] == 3:
                            saMask[0:saSidePix,int(pupilSidePix/2-saSidePix/2):int(pupilSidePix/2+saSidePix/2)] = 1
                            saMask *= cpuArray(self.fao.ao.tel.pupil)
                        else:
                            saMask[int(pupilSidePix/2-saSidePix/2):int(pupilSidePix/2+saSidePix/2),\
                                   int(pupilSidePix/2-saSidePix/2):int(pupilSidePix/2+saSidePix/2)] = 1
                        if gpuMastsel:
                            maskFocus.sampling = congrid(cp.asarray(saMask), [self.sx, self.sx])
                        else:
                            maskFocus.sampling = congrid(saMask, [self.sx, self.sx])
                        maskFocus.sampling = zeroPad(maskFocus.sampling, (self.N-self.sx)//2)
                        if nMaskFocus > 1:
                            self.maskFocus.append(maskFocus)
                        else:
                            self.maskFocus = maskFocus

                # -----------------------------------------------------------------
                ## PSF for NGS directions

                if self.verbose:
                    print('******** Focus Sensor PSF - NGS directions (1 sub-aperture)')

                Focus_SR, psdArray, psfLE_Focus, Focus_FWHM_mas = self.psdSetToPsfSet(self.N,
                                                                           self.freq_range,
                                                                           self.dk,
                                                                           self.maskFocus,
                                                                           PSD_Focus,
                                                                           self.Focus_wvl,
                                                                           self.psInMas[0],
                                                                           self.nPixPSF,
                                                                           scaleFactor=(2*np.pi*1e-9/self.Focus_wvl)**2,
                                                                           oversampling=self.oversampling)

                # -----------------------------------------------------------------
                ## Merit functions
                self.Focus_SR_field         = Focus_SR
                self.Focus_FWHM_mas_field   = Focus_FWHM_mas
                self.Focus_EE_field         = []
                idx = 0
                for img in psfLE_NGS:
                    ee_,rr_ = getEncircledEnergy(img.sampling, pixelscale=self.psInMas[0], center=(self.fao.ao.cam.fovInPix/2,self.fao.ao.cam.fovInPix/2), nargout=2)
                    ee_ *= 1/np.max(ee_)
                    ee_at_radius_fn = interp1d(rr_, ee_, kind='cubic', bounds_error=False)
                    # max is used to compute EE on at least a radius of one pixel
                    ee_Focus = ee_at_radius_fn(max([Focus_FWHM_mas[idx],self.Focus_psInMas]))
                    self.Focus_EE_field.append(ee_Focus)
                    if self.verbose:
                        print('EE (focus sensor)     :', "%.5f" % ee_Focus)
                    idx += 1
            else:
                if self.verbose:
                    print('Focus sensor is not set: using LO PSFs.')
                self.Focus_SR_field         = self.NGS_SR_field
                self.Focus_FWHM_mas_field   = self.NGS_FWHM_mas_field
                self.Focus_EE_field         = self.NGS_EE_field


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
            self.sr.append( np.exp( -4*np.pi**2 * ( self.penalty[-1]**2 )/(self.wvl*1e9)**2) )
            scale = (np.pi/(180*3600*1000) * self.TelescopeDiameter / (4*1e-9))
            fwhms_lo = 2.355 * self.LO_res/scale / np.sqrt(2)
            self.fwhm.append(np.sqrt(fwhms_lo**2 + np.asarray(self.pointings_FWHM_mas)**2))
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
            for img in self.results:
                self.sr.append(getStrehl(img.sampling, self.fao.ao.tel.pupil, self.fao.freq.sampRef, method='otf'))
                self.fwhm.append(getFWHM(img.sampling, self.psInMas[0], method='contour', nargout=1))
                if self.ensquaredEnergy:
                    ee_ = cpuArray(getEnsquaredEnergy(self.cubeResultsArray[i,:,:]))
                    rr_ = np.arange(1, ee_.shape[0]*2, 2) * self.psInMas[0] * 0.5
                else:
                    ee_,rr_ = getEncircledEnergy(img.sampling, pixelscale=self.psInMas[0],
                                                 center=(self.fao.ao.cam.fovInPix/2,self.fao.ao.cam.fovInPix/2), nargout=2)
                ee_at_radius_fn = interp1d(rr_, ee_, kind='cubic', bounds_error=False)
                self.ee.append( cpuArray(ee_at_radius_fn(self.eeRadiusInMas)).item() )


    def doOverallSimulation(self, astIndex=None):

        if self.LOisOn:        
            self.configLO(astIndex)
        
        self.results = []

        # ------------------------------------------------------------
        # **** Start of calculations only performed on first call ****
        # ------------------------------------------------------------

        if astIndex is None or self.firstSimCall:
            # ------------------------------------------------------------------------
            ## HO Part with P3 PSDs

            if self.verbose:
                print('******** HO PSD science and NGSs directions')

            self.fao = fourierModel( self.fullPathFilename, calcPSF=False, verbose=self.verbose
                               , display=False, getPSDatNGSpositions=self.LOisOn
                               , computeFocalAnisoCov=False, TiltFilter=self.LOisOn
                               , getErrorBreakDown=self.getHoErrorBreakDown, doComputations=False)

            if 'sensor_LO' in self.my_data_map.keys():
                self.fao.my_data_map['sensor_LO']['NumberPhotons'] = self.my_data_map['sensor_LO']['NumberPhotons']
                self.fao.ao.my_data_map['sensor_LO']['NumberPhotons'] = self.my_data_map['sensor_LO']['NumberPhotons']
            if 'sources_LO' in self.my_data_map.keys():
                self.fao.my_data_map['sources_LO'] = self.my_data_map['sources_LO']
                self.fao.ao.my_data_map['sources_LO'] = self.my_data_map['sources_LO']
                self.fao.ao.configLOsensor()
                self.fao.ao.configLO()
                self.fao.ao.configLO_SC()

            self.fao.initComputations()

            # High-order PSD caculations at the science directions and NGSs directions
            self.PSD           = self.fao.PSD # in nm^2
            self.PSD           = self.PSD.transpose()
            self.N             = self.PSD[0].shape[0]
            self.nPointings    = self.pointings.shape[1]
            self.nPixPSF       = int(self.fao.freq.nOtf /self.fao.freq.kRef_)
            if isinstance(self.fao.freq.k_, list):
                self.oversampling = self.fao.freq.k_[0]
            else:
                self.oversampling  = self.fao.freq.k_
            self.freq_range    = self.fao.ao.cam.fovInPix*self.fao.freq.PSDstep*self.oversampling
            self.pitch         = 1/self.freq_range
            self.grid_diameter = self.pitch*self.N
            self.sx            = int(2*np.round(self.tel_radius/self.pitch))
            self.dk            = 1e9*self.fao.freq.kcMax_/self.fao.freq.resAO
            # Define the pupil shape
            self.mask = Field(self.wvl, self.N, self.grid_diameter)
            self.psInMas = cpuArray(self.fao.freq.psInMas)
            self.mask.sampling = congrid(arrayP3toMastsel(self.fao.ao.tel.pupil), [self.sx, self.sx])
            self.mask.sampling = zeroPad(self.mask.sampling, (self.N-self.sx)//2)

            if self.verbose:
                print('fao.samp:', self.fao.freq.samp)
                print('fao.PSD.shape:', self.fao.PSD.shape)
                print('fao.freq.psInMas:', self.psInMas)

            # ------------------------------------------------------------------------
            ## Update the instrumental OTF if static WFE is present:
            ##     [telescope] PathStaticOn key in params
            if (self.fao.ao.tel.opdMap_on is not None):
                self.otfStatic = []
                # following lines come from SF2PSF method of P3/p3/aoSystem/FourierUtils.py
                for jWvl in range(self.fao.freq.nWvl):
                    otfStatic, otfDL, phaseMap = \
                    getStaticOTF(self.fao.ao.tel,int(self.fao.freq.nOtf),self.fao.freq.samp[jWvl],
                                 self.fao.freq.wvl[jWvl],spatialFilter=1)
                    self.otfStatic.append(arrayP3toMastsel(otfStatic))

            # ----------------------------------------------------------------------------
            ## optional LO part
            if self.LOisOn:
                if self.verbose:
                    print('******** LO PART')
                    
                # ------------------------------------------------------------------------
                # --- NGS PSDs, PSFs and merit functions on PSFs
                self.ngsPSF()
                
                # ------------------------------------------------------------------------
                # --- initialize MASTSEL MavisLO object
                self.mLO = MavisLO(self.path, self.parametersFile, verbose=self.verbose)

        # ------------------------------------------------------------
        # ***** End of calculations only performed on first call *****
        # ------------------------------------------------------------

        # ----------------------------------------------------------------------------
        # optional LO part - for every call, not just the first (self.firstSimCall)
        if self.LOisOn:
            # ------------------------------------------------------------------------
            ## total covariance matrix Ctot
            if astIndex is None:
                self.Ctot          = self.mLO.computeTotalResidualMatrix(np.array(self.cartSciencePointingCoords),
                                                                         self.cartNGSCoords_field, self.NGS_fluxes_field,
                                                                         self.LO_freqs_field,
                                                                         self.NGS_SR_field, self.NGS_EE_field, self.NGS_FWHM_mas_field, doAll=True)

                # --------------------------------------------------------------------
                # --- optional total focus covariance matrix Ctot
                if self.addFocusError:

                    # compute focus error
                    self.CtotFocus = self.mLO.computeFocusTotalResidualMatrix(self.cartNGSCoords_field, self.Focus_fluxes_field,
                                                                         self.Focus_freqs_field, self.Focus_SR_field,
                                                                         self.Focus_EE_field, self.Focus_FWHM_mas_field)

                    self.GF_res = np.sqrt(self.CtotFocus[0])

                    # add focus error to PSD using P3 FocusFilter
                    FocusFilter = self.fao.FocusFilter()
                    FocusFilter *= 1/FocusFilter.sum()
                    for PSDho in self.PSD:
                        PSDho += self.GF_res**2 * FocusFilter
                    self.GFinPSD = True
                # ---------------------------------------------------------------------
            else:                  
                if self.firstSimCall:
                    self.mLO.computeTotalResidualMatrix(np.array(self.cartSciencePointingCoords),
                                                        self.cartNGSCoords_field, self.NGS_fluxes_field,
                                                        self.LO_freqs_field, self.NGS_SR_field,
                                                        self.NGS_EE_field, self.NGS_FWHM_mas_field, doAll=False)
                    if self.addFocusError:
                        self.mLO.computeFocusTotalResidualMatrix(self.cartNGSCoords_field, self.Focus_fluxes_field,
                                                                 self.Focus_freqs_field, self.Focus_SR_field,
                                                                 self.Focus_EE_field, self.Focus_FWHM_mas_field)

                # discard guide stars with flux less than 1 photon per frame per subaperture
                if np.min(self.NGS_fluxes_asterism) < 1 and np.max(self.NGS_fluxes_asterism) > 1:
                    valid_indices = np.where(np.array(self.NGS_fluxes_asterism) > 1)[0]
                    self.NGS_fluxes_asterism = [elem for i, elem in enumerate(self.NGS_fluxes_asterism) if i in valid_indices]
                    self.Focus_fluxes_asterism = [elem for i, elem in enumerate(self.Focus_fluxes_asterism) if i in valid_indices]
                    self.cartNGSCoords_asterism = [elem for i, elem in enumerate(self.cartNGSCoords_asterism) if i in valid_indices]
                    self.currentAsterismIndices = [elem for i, elem in enumerate(self.currentAsterismIndices) if i in valid_indices]

                self.Ctot  = self.mLO.computeTotalResidualMatrixI(self.currentAsterismIndices,
                                                                  np.array(self.cartSciencePointingCoords),
                                                                  np.array(self.cartNGSCoords_asterism),
                                                                  self.NGS_fluxes_asterism)
                
                # --------------------------------------------------------------------
                # --- optional total focus covariance matrix Ctot
                if self.addFocusError:
                    self.CtotFocus = self.mLO.computeFocusTotalResidualMatrixI(self.currentAsterismIndices,
                                                                         np.array(self.cartNGSCoords_asterism),
                                                                         self.Focus_fluxes_asterism)
                    self.GF_res = np.sqrt(self.CtotFocus[0])
                    self.GFinPSD = False
                # ---------------------------------------------------------------------

            # ------------------------------------------------------------------------
            ## computation of the LO error (this changes for each asterism)
            self.LO_res = np.sqrt(np.trace(self.Ctot,axis1=1,axis2=2))
        # ------------------------------------------------------------------------

        # ------------------------------------------------------------------------
        # final PSF computation
        self.finalPSF(astIndex)

        # ------------------------------------------------------------------------
        # plots
        if self.doPlot:
            if self.LOisOn and self.doConvolve:
                tiledDisplay(self.results)
                plotEllipses(self.cartSciencePointingCoords, self.cov_ellipses, 0.4)
            else:
                self.results[0].standardPlot(True)

        # set first call attribute to False
        self.firstSimCall = False

        # ------------------------------------------------------------------------
        # final results
        if astIndex is None:
            self.computeOL_PSD()
            self.computeDL_PSD()
            self.cubeResults = []
            for img in self.results:
                self.cubeResults.append(img.sampling)
            self.cubeResultsArray = np.array(self.cubeResults)

            if self.verbose:
                print('HO_res [nm]:',self.HO_res)
                if self.LOisOn:
                    print('LO_res [nm]:',self.LO_res)
                if hasattr(self,'GF_res'):
                    print('GF_res [nm]:',self.GF_res)
