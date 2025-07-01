import p3.aoSystem
import p3.aoSystem.fourierModel
import p3.aoSystem.FourierUtils

from p3.aoSystem.fourierModel import *
from p3.aoSystem.FourierUtils import *
from mastsel import *

from .tiptopUtils import *
from ._version import __version__

from matplotlib import cm
import matplotlib as mpl
norm = mpl.colors.Normalize(vmin=0, vmax=1)
rc("text", usetex=False)

import json

rad2mas = 3600 * 180 * 1000 / np.pi

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
                          doPlot=False, addSrAndFwhm=True,
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
            self.wvlMax = max(wvl_temp)
            self.wvl = wvl_temp     # lambda
            self.nWvl = len(wvl_temp)
        else:
            self.wvlMax = wvl_temp
            self.wvl = [wvl_temp]     # lambda
            self.nWvl = 1
        self.zenithSrc  = self.my_data_map['sources_science']['Zenith']
        self.azimuthSrc = self.my_data_map['sources_science']['Azimuth']
        self.pointings = polarToCartesian(np.array( [self.zenithSrc, self.azimuthSrc]))
        self.xxSciencePointigs         = self.pointings[0,:]
        self.yySciencePointigs         = self.pointings[1,:]
        self.psInMas = self.my_data_map['sensor_science']['PixelScale']
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

        self.LO_zen_field       = self.my_data_map['sources_LO']['Zenith']
        self.LO_az_field        = self.my_data_map['sources_LO']['Azimuth']
        self.LO_fluxes_field    = self.my_data_map['sensor_LO']['NumberPhotons']
        self.LO_psInMas         = self.my_data_map['sensor_LO']['PixelScale']
        # if self.LO_psInMas is a scalar makes a list on n elements
        if not isinstance(self.LO_psInMas, list):
            self.LO_psInMas = [self.LO_psInMas] * len(self.LO_zen_field)
        self.LO_freqs_field     = self.my_data_map['RTC']['SensorFrameRate_LO']
        if not isinstance(self.LO_freqs_field, list):
            self.LO_freqs_field = [self.LO_freqs_field] * len(self.LO_zen_field)
        if self.check_config_key('sensor_LO','addAliasError'):
            self.addLoAlias     = self.my_data_map['sensor_LO']['addAliasError']
        else:
            self.addLoAlias     = False

        if self.check_section_key('sensor_Focus'):
            self.Focus_fluxes4s_field   = self.my_data_map['sensor_Focus']['NumberPhotons']
            self.Focus_psInMas          = self.my_data_map['sensor_Focus']['PixelScale']
            # if self.Focus_psInMas is a scalar makes a list on n elements
            if not isinstance(self.Focus_psInMas, list):
                self.Focus_psInMas = [self.Focus_psInMas] * len(self.LO_zen_field)
            if self.check_section_key('sources_Focus'):
                Focus_wvl_temp          = self.my_data_map['sources_Focus']['Wavelength']
            else:
                Focus_wvl_temp          = self.my_data_map['sources_LO']['Wavelength']
            if isinstance(Focus_wvl_temp, list):
                self.Focus_wvl          = Focus_wvl_temp[0]  # lambda
            else:
                self.Focus_wvl          = Focus_wvl_temp     # lambda
        else:
            self.Focus_fluxes4s_field   = self.LO_fluxes_field
            self.Focus_psInMas          = self.LO_psInMas
            self.Focus_wvl              = self.LO_wvl
        if self.check_config_key('RTC','SensorFrameRate_Focus'):
            self.Focus_freqs_field      = self.my_data_map['RTC']['SensorFrameRate_Focus']
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
        polarNGSCoords = np.asarray(polarNGSCoordsList)
        self.nNaturalGS_field = len(self.LO_zen_field)
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

    def computePSF1D(self):
        psf1d = []
        psf1d_radius = None
        psf1d_radius_list_list = []
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
                rr, radialprofile, ee = radial_profile(psf, ext=0, pixelscale=self.psInMas, ee=True,
                                                       center=center, stddev=False, binsize=None, maxradius=self.psInMas*psfRadius,
                                                       normalize='total', pa_range=None, slice=0, nargout=2, verbose=self.verbose)
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
        jsondict['infos'] = self.my_data_map
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
                hdr1['RESH'+str(i).zfill(4)] =  np.round(cpuArray(self.HO_res[i]),3)
        if hasattr(self,'LO_res'):
            hdr1['RESL'] = "Low Order residual in nm RMS"
            for i in range(self.LO_res.shape[0]):
                hdr1['RESL'+str(i).zfill(4)] = np.round(cpuArray(self.LO_res[i]),3)
        if hasattr(self,'GF_res'):
            hdr1['RESF'] = "Global Focus residual in nm RMS (included in PSD)"
            hdr1['RESF0000'] = np.round(cpuArray(self.GF_res),3)
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
                    hdr1[fTxt+str(j).zfill(Nfill)+wTxt] = np.round(fwhm_temp,3)
                for j in range(cubeResultsArray.shape[0]):
                    if self.ensquaredEnergy:
                        ee = cpuArray(getEnsquaredEnergy(cubeResultsArray[j,:,:]))
                        rr = np.arange(1, ee.shape[0]*2, 2) * self.psInMas * 0.5
                    else:
                        ee,rr = getEncircledEnergy(cubeResultsArray[j,:,:], pixelscale=self.psInMas,
                                                   center=(self.nPixPSF/2,self.nPixPSF/2), nargout=2)
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

        hdul1.writeto( os.path.join(self.outputDir, self.outputFile + '.fits'), overwrite=True)
        if self.verbose:
            print("Output cube shape:", self.cubeResultsArray.shape)

    def computeOL_PSD(self):
        # OPEN-LOOP PSD
        k = np.sqrt(self.fao.freq.k2_)
        pf = FourierUtils.pistonFilter(2*self.tel_radius,k)
        spectrum = arrayP3toMastsel(self.fao.ao.atm.spectrum(k) * pf)
        psdOL = Field(self.wvlRef, self.N, self.freq_range, 'rad')
        psdOL.sampling = spectrum * (self.dk*self.wvlRef/np.pi)**2 # the PSD must be provided in m^2.m^2
        padPSD = self.nWvl > 1
        mask = arrayP3toMastsel(self.fao.ao.tel.pupil)
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


    def finalConvolution(self):
        self.cov_ellipses = self.mLO.ellipsesFromCovMats(self.Ctot)
        if self.verbose:
            for n in range(self.cov_ellipses.shape[0]):
                print('cov_ellipses #',n,': ', self.cov_ellipses[n,:], ' (unit: rad, mas, mas)')
            print('******** FINAL CONVOLUTION')
        # CONVULUTION KERNELS
        resSpecList = []
        resSpecListJ = []
        for ellp in self.cov_ellipses:
            resSpecList.append(residualToSpectrum(ellp, self.wvlRef, self.nPixPSF, 1/(self.nPixPSF * self.psInMas)))
            if self.jitter_FWHM is not None:
                if isinstance(self.jitter_FWHM, list):
                    ellpJ = [self.jitter_FWHM[2], sigma_from_FWHM(self.jitter_FWHM[0]), sigma_from_FWHM(self.jitter_FWHM[1])]
                else:
                    ellpJ = [0, sigma_from_FWHM(self.jitter_FWHM), sigma_from_FWHM(self.jitter_FWHM)]
                resSpecListJ.append(residualToSpectrum(ellpJ, self.wvlRef, self.nPixPSF, 1/(self.nPixPSF * self.psInMas)))
            else:
                resSpecListJ.append(0)
        # FINAl CONVOLUTION
        for i in range(self.nWvl):
            if self.nWvl>1:
                psfList = self.psfLongExpPointingsArr[i]
            else:
                psfList = self.psfLongExpPointingsArr
            resultList = []
            for psfLongExp, resSpec, resSpecJ in zip(psfList, resSpecList, resSpecListJ):
                temp = convolve(psfLongExp, resSpec)
                if self.jitter_FWHM is not None:
                    temp = convolve(temp, resSpecJ)
                resultList.append(temp)
            if self.nWvl>1:
                self.results.append(resultList)
            else:
                self.results = resultList

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
                if self.nWvl>1:
                    self.pointings_FWHM_mas.append(fwhmList)
                else:
                    self.pointings_FWHM_mas = fwhmList

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
                for i in range(self.nWvl):
                    if self.nWvl>1:
                        psfList = self.psfLongExpPointingsArr[i]
                    else:
                        psfList = self.psfLongExpPointingsArr
                    resultList = []
                    for psfLongExp in psfList:
                        resultList.append(psfLongExp)
                    if self.nWvl>1:
                        self.results.append(resultList)
                    else:
                        self.results = resultList
        else:
            if self.jitter_FWHM is not None:
                if isinstance(self.jitter_FWHM, list):
                    ellp = [self.jitter_FWHM[2], sigma_from_FWHM(self.jitter_FWHM[0]), sigma_from_FWHM(self.jitter_FWHM[1])]
                else:
                    ellp = [0, sigma_from_FWHM(self.jitter_FWHM), sigma_from_FWHM(self.jitter_FWHM)]
                resSpecJ = residualToSpectrum(ellp, self.wvlRef, self.nPixPSF, 1/(self.nPixPSF * self.psInMas))
            for i in range(self.nWvl):
                if self.nWvl>1:
                    psfList = self.psfLongExpPointingsArr[i]
                else:
                    psfList = self.psfLongExpPointingsArr
                resultList = []
                for psfLongExp in psfList:
                    if self.jitter_FWHM is not None:
                        resultList.append(convolve(psfLongExp,resSpecJ))
                    else:
                        resultList.append(psfLongExp)
                if self.nWvl>1:
                    self.results.append(resultList)
                else:
                    self.results = resultList

    def ngsPSF(self):
        # pixel size for LO
        LO_PSFsInMas = self.psInMas*self.LO_wvl/self.wvlMax

        # skip reshape in psdSetToPsfSet to get a high sampling PSF if original sampling is low
        if LO_PSFsInMas/np.min(self.LO_psInMas) > 1 and self.overSamp > 1:
            skip_reshape = True
            LO_PSFsInMas /= self.overSamp
            nPixPSFLO = int(self.overSamp * self.nPixPSF)
        else:
            skip_reshape = False
            nPixPSFLO = self.nPixPSF

        # -----------------------------------------------------------------
        # PSD and sub-aperture mask for NGS directions
        psdNGS = arrayP3toMastsel(self.PSD[-self.nNaturalGS_field:])
        k  = np.sqrt(self.fao.freq.k2_)

        # Define the LO sub-aperture shape
        nSA = self.my_data_map['sensor_LO']['NumberLenslets']
        maskLO = maskSA(nSA, self.nNaturalGS_field, arrayP3toMastsel(self.fao.ao.tel.pupil))

        for i in range(self.nNaturalGS_field):
            # This is needed when the AsterismSelection has many more stars than the length of nSA.
            if len(nSA) == self.nNaturalGS_field:
                nSAi = nSA[i]
                len_nSA = len(nSA)
            else:
                nSAi = nSA[0]
                len_nSA = 1
            if nSAi != 1:
                # piston filter for the sub-aperture size
                pf = FourierUtils.pistonFilter(2*self.tel_radius/nSAi,k)
                psdNGS[i] = psdNGS[i] * pf

        # -----------------------------------------------------------------
        # PSF for NGS directions

        if self.verbose:
            print('******** LO PSF - NGS directions (1 sub-aperture)')
        psfLE_NGS = psdSetToPsfSet(psdNGS, maskLO,
                                   self.LO_wvl, self.N, self.sx, self.grid_diameter,
                                   self.freq_range, self.dk, nPixPSFLO,
                                   self.wvlMax, self.overSamp,
                                   opdMap=self.opdMap, skip_reshape=skip_reshape)

        # -----------------------------------------------------------------
        # Merit functions
        self.NGS_SR_field           = []
        self.NGS_FWHM_mas_field     = []
        self.NGS_EE_field           = []
        idx = 0
        for img in psfLE_NGS:
            # Get SR, FWHM in mas and EE at the NGSs positions at the sensing wavelength
            s1 = cpuArray(psdNGS[idx]).sum()
            SR = np.exp(-s1*(2*np.pi*1e-9/self.LO_wvl)**2) # Strehl-ratio at the sensing wavelength
            self.NGS_SR_field.append(SR)
            fwhmX,fwhmY = getFWHM(img.sampling, LO_PSFsInMas, method='contour', nargout=2)
            FWHM = np.sqrt(fwhmX*fwhmY) #average over major and minor axes
            self.NGS_FWHM_mas_field.append(FWHM)
            if 2*FWHM >= nPixPSFLO*LO_PSFsInMas:
                ee_NGS = 1
            else:
                ee_,rr_ = getEncircledEnergy(img.sampling, pixelscale=LO_PSFsInMas,
                                             center=(nPixPSFLO/2,nPixPSFLO/2), nargout=2)
                ee_ *= 1/np.max(ee_)
                ee_at_radius_fn = interp1d(rr_, ee_, kind='cubic', bounds_error=False)
                # max is used to compute EE on at least a radius of one pixel
                if isinstance(self.LO_psInMas,list):
                    LO_psInMas_i = self.LO_psInMas[idx]
                else:
                    LO_psInMas_i = self.LO_psInMas
                ee_NGS = ee_at_radius_fn(max([FWHM,LO_psInMas_i]))
            self.NGS_EE_field.append(ee_NGS)
            if self.verbose:
                print('SR(@',int(self.LO_wvl*1e9),'nm)        :', "%.5f" % SR)
                print('FWHM(@',int(self.LO_wvl*1e9),'nm) [mas]:', "%.3f" % FWHM)
                print('EE                   :', "%.5f" % ee_NGS)
            idx += 1

        if self.addLoAlias:
            if self.verbose:
                print('Adding aliasing error on LO!')
            # DIFFRACTION LIMITED PSD and PSF
            if isinstance(maskLO, list):
                self.NGS_DL_FWHM_mas = []
            else:
                self.NGS_DL_FWHM_mas = None
            for i in range(len_nSA):
                if isinstance(maskLO, list):
                    maskI = maskLO[i]
                else:
                    maskI = maskLO
                psdDL = Field(self.LO_wvl, self.N, self.freq_range, 'rad')
                maskField = Field(self.LO_wvl, self.N, self.grid_diameter)
                maskField.sampling = congrid(maskI, [self.sx, self.sx])
                maskField.sampling = zeroPad(maskField.sampling, (self.N-self.sx)//2)
                psfNgsDL = longExposurePsf(maskField, psdDL)
                fwhmX,fwhmY  = getFWHM( psfNgsDL.sampling, LO_PSFsInMas, method='contour', nargout=2)
                if self.NGS_DL_FWHM_mas is None:
                    self.NGS_DL_FWHM_mas = np.sqrt(fwhmX*fwhmY)
                else:
                    self.NGS_DL_FWHM_mas.append(np.sqrt(fwhmX*fwhmY))
        else:
            self.NGS_DL_FWHM_mas = None

        # -----------------------------------------------------------------
        # optional Focus error
        if self.addFocusError:
            # pixel size for Focus
            Focus_PSFsInMas = self.psInMas*self.Focus_wvl/self.wvlMax

            # skip reshape in psdSetToPsfSet to get a high sampling PSF if original sampling is low
            if Focus_PSFsInMas/np.min(self.Focus_psInMas) > 1 and self.overSamp > 1:
                skip_reshape = True
                Focus_PSFsInMas /= self.overSamp
                nPixPSFFocus = int(self.overSamp * self.nPixPSF)
            else:
                skip_reshape = False
                nPixPSFFocus = self.nPixPSF

            if 'sensor_Focus' in self.my_data_map.keys():
                if self.verbose:
                    print('Focus sensor is set: computing new PSFs.')
                # -----------------------------------------------------------------
                ## PSD and sub-aperture mask for NGS directions

                # Define the LO sub-aperture shape
                nSAfocus = self.my_data_map['sensor_Focus']['NumberLenslets']
                psdFocus = arrayP3toMastsel(self.PSD[-self.nNaturalGS_field:])
                maskFocus = maskSA(nSAfocus, self.nNaturalGS_field, arrayP3toMastsel(self.fao.ao.tel.pupil))

                for i in range(self.nNaturalGS_field):
                    # This is needed when the AsterismSelection has many more stars than the length of nSAfocus.
                    if len(nSAfocus) == self.nNaturalGS_field:
                        nSAfocusI = nSAfocus[i]
                    else:
                        nSAfocusI = nSAfocus[0]
                    if nSAfocusI != 1:
                        # --- piston filter for the sub-aperture size
                        pf = FourierUtils.pistonFilter(2*self.tel_radius/nSAfocusI,k)
                        psdFocus[i] = psdFocus[i] * pf


                # -----------------------------------------------------------------
                ## PSF for NGS directions

                if self.verbose:
                    print('******** Focus Sensor PSF - NGS directions (1 sub-aperture)')
                psfLE_Focus = psdSetToPsfSet(psdFocus, maskFocus,
                                             self.Focus_wvl, self.N, self.sx, self.grid_diameter,
                                             self.freq_range, self.dk, nPixPSFFocus,
                                             self.wvlMax, self.overSamp,
                                             opdMap=self.opdMap, skip_reshape=skip_reshape)

                # -----------------------------------------------------------------
                ## Merit functions
                self.Focus_SR_field         = []
                self.Focus_FWHM_mas_field   = []
                self.Focus_EE_field         = []
                idx = 0
                for img in psfLE_Focus:
                    # Get SR, FWHM in mas and EE at the NGSs positions at the sensing wavelength
                    s1 = cpuArray(psdFocus[idx]).sum()
                    SR = np.exp(-s1*(2*np.pi*1e-9/self.Focus_wvl)**2) # Strehl-ratio at the sensing wavelength
                    self.Focus_SR_field.append(SR)
                    fwhmX,fwhmY = getFWHM(img.sampling, Focus_PSFsInMas, method='contour', nargout=2)
                    FWHM = np.sqrt(fwhmX*fwhmY) #average over major and minor axes
                    self.Focus_FWHM_mas_field.append(FWHM)
                    if 2*FWHM >= nPixPSFFocus*Focus_PSFsInMas:
                        ee_Focus = 1
                    else:
                        ee_,rr_ = getEncircledEnergy(img.sampling, pixelscale=Focus_PSFsInMas,
                                                     center=(nPixPSFFocus/2,nPixPSFFocus/2), nargout=2)
                        ee_ *= 1/np.max(ee_)
                        ee_at_radius_fn = interp1d(rr_, ee_, kind='cubic', bounds_error=False)
                        # max is used to compute EE on at least a radius of one pixel
                        if isinstance(self.Focus_psInMas,list):
                            Focus_psInMas_i = self.Focus_psInMas[idx]
                        else:
                            Focus_psInMas_i = self.Focus_psInMas
                        ee_Focus = ee_at_radius_fn(max([FWHM,Focus_psInMas_i]))
                    self.Focus_EE_field.append(ee_Focus)
                    if self.verbose:
                        print('SR(@',int(self.Focus_wvl*1e9),'nm)        :', "%.5f" % SR)
                        print('FWHM(@',int(self.Focus_wvl*1e9),'nm) [mas]:', "%.3f" % FWHM)
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
            self.sr.append( np.exp( -4*np.pi**2 * ( self.penalty[-1]**2 )/(self.wvlRef*1e9)**2) )
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

            for i in range(self.nWvl):
                if self.nWvl>1:
                    results = self.results[i]
                else:
                    results = self.results
                samp = self.wvl[i] * rad2mas / (self.psInMas*2*self.tel_radius)
                sr = []
                fwhm = []
                ee = []
                for img in results:
                    sr.append(getStrehl(img.sampling, self.fao.ao.tel.pupil, samp, method='max', psfInOnePix=True))
                    fwhm.append(getFWHM(img.sampling, self.psInMas, method='contour', nargout=1))
                    if self.ensquaredEnergy:
                        ee_ = cpuArray(getEnsquaredEnergy(img.sampling))
                        rr_ = np.arange(1, ee_.shape[0]*2, 2) * self.psInMas * 0.5
                    else:
                        ee_,rr_ = getEncircledEnergy(img.sampling, pixelscale=self.psInMas,
                                                     center=(self.nPixPSF/2,self.nPixPSF/2), nargout=2)
                    ee_at_radius_fn = interp1d(rr_, ee_, kind='cubic', bounds_error=False)
                    ee.append( cpuArray(ee_at_radius_fn(self.eeRadiusInMas)).item() )
                if self.nWvl>1:
                    self.sr.append(sr)
                    self.fwhm.append(fwhm)
                    self.ee.append(ee)
                else:
                    self.sr = sr
                    self.fwhm = fwhm
                    self.ee = ee


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
                               , getErrorBreakDown=self.getHoErrorBreakDown, doComputations=False
                               , psdExpansion=True, reduce_memory=True)

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
            self.nPixPSF       = self.my_data_map['sensor_science']['FieldOfView']
            self.overSamp      = int(self.fao.freq.kRef_)
            self.PSDstep       = self.fao.freq.PSDstep
            self.freq_range    = self.N*self.PSDstep
            self.grid_diameter = 1/self.PSDstep
            self.sx            = int(2*np.round(self.tel_radius*self.freq_range))
            # dk is the same as in p3.aoSystem.powerSpectrumDensity except that it is multiplied by 1e9 instead of 2.
            self.dk            = 1e9*self.fao.freq.kcMax_/self.fao.freq.resAO
            # wvlRef from P3 is required to scale correctly the OL PSD from rad to m
            self.wvlRef        = self.fao.freq.wvlRef
            # Define the pupil shape
            self.mask = Field(self.wvlRef, self.N, self.grid_diameter)
            self.mask.sampling = congrid(arrayP3toMastsel(self.fao.ao.tel.pupil), [self.sx, self.sx])
            self.mask.sampling = zeroPad(self.mask.sampling, (self.N-self.sx)//2)
            # error messages for wrong pixel size
            if self.psInMas != cpuArray(self.fao.freq.psInMas[0]):
                raise ValueError("sensor_science.PixelScale, '{}', is different from self.fao.freq.psInMas,'{}'"
                         .format(self.psInMas,cpuArray(self.fao.freq.psInMas)))

            if self.fao.ao.tel.opdMap_on is not None:
                self.opdMap = arrayP3toMastsel(self.fao.ao.tel.opdMap_on)
            else:
                self.opdMap = None

            if self.verbose:
                print('PSD step:', self.PSDstep)
                print('PSD freq range:', self.freq_range)
                print('PSD shape:', self.PSD.shape)
                print('oversampling:', self.overSamp)
                print('sensor_science.PixelScale:', self.psInMas)

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
                                                                         self.NGS_SR_field, self.NGS_EE_field, self.NGS_FWHM_mas_field,
                                                                         aNGS_FWHM_DL_mas = self.NGS_DL_FWHM_mas, doAll=True)

                # --------------------------------------------------------------------
                # --- optional total focus covariance matrix Ctot
                if self.addFocusError:

                    # compute focus error
                    self.CtotFocus = self.mLO.computeFocusTotalResidualMatrix(self.cartNGSCoords_field, self.Focus_fluxes_field,
                                                                         self.Focus_freqs_field, self.Focus_SR_field,
                                                                         self.Focus_EE_field, self.Focus_FWHM_mas_field)

                    self.GF_res = np.sqrt(max(self.CtotFocus[0], 0))
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
                                                        self.NGS_EE_field, self.NGS_FWHM_mas_field,
                                                        aNGS_FWHM_DL_mas = self.NGS_DL_FWHM_mas, doAll=False)
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
                    self.GF_res = np.sqrt(max(self.CtotFocus[0], 0))
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
            if self.nWvl>1:
                results = self.results[0]
            else:
                results = self.results
            if self.LOisOn and self.doConvolve:
                tiledDisplay(results)
                plotEllipses(self.cartSciencePointingCoords, self.cov_ellipses, 0.4)
            else:
                results[0].standardPlot(True)

        # set first call attribute to False
        self.firstSimCall = False

        # ------------------------------------------------------------------------
        # final results
        if astIndex is None:
            self.computeOL_PSD()
            self.computeDL_PSD()
            self.cubeResults = []
            cubeResultsArray = []
            for i in range(self.nWvl):
                if self.nWvl>1:
                    results = self.results[i]
                else:
                    results = self.results
                cubeResults = []
                for img in results:
                    cubeResults.append(cpuArray(img.sampling))
                if self.nWvl>1:
                    self.cubeResults.append(cubeResults)
                    cubeResultsArray.append(np.array(cubeResults))
                else:
                    self.cubeResults = cubeResults
                    cubeResultsArray = cubeResults
                self.cubeResultsArray = np.array(cubeResultsArray)
            self.computePSF1D()
            if self.verbose:
                print('HO_res [nm]:',self.HO_res)
                if self.LOisOn:
                    print('LO_res [nm]:',self.LO_res)
                if hasattr(self,'GF_res'):
                    print('GF_res [nm]:',self.GF_res)
