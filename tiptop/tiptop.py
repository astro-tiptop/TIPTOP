import os
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import rc
from p3.aoSystem.fourierModel import *
from p3.aoSystem.FourierUtils import *
from configparser import ConfigParser
import yaml

from mastsel import *
from mastsel import gpuEnabled as gpuMastsel
from p3.aoSystem import gpuEnabled as gpuP3

from datetime import datetime

rc("text", usetex=False)

def arrayP3toMastsel(v):
    if (gpuP3 and gpuMastsel) or (not gpuP3 and not gpuMastsel):
        return v
    elif not gpuP3 and gpuMastsel:
        return cp.asarray(v)
    elif gpuP3 and not gpuMastsel: 
        return v.get()
    
def cpuArray(v):
    if isinstance(v,np.ndarray):
        return v
    else:
        return v.get()

MAX_VALUE_CHARS = 80
APPEND_TOKEN = '&&&'

def add_hdr_keyword(hdr, key_primary, key_secondary, val, iii=None, jjj=None):
    val_string = str(val)
    key = 'HIERARCH '+ key_primary +' '+ key_secondary
    if not iii is None:
        key += ' '+str(iii)    
    if not jjj is None:
        key +=  ' '+str(jjj)    
    margin = 4
    key = key
    current_val_string = val_string
    if len(key) + margin > MAX_VALUE_CHARS:
        print("Error, keywork is not acceptable due to string length.")
        return
    while not len(key) + 1 + len(current_val_string) + margin < MAX_VALUE_CHARS:
        max_char_index = MAX_VALUE_CHARS-len(key)-1-len(APPEND_TOKEN)-margin
        hdr[key+'+'] = current_val_string[:max_char_index]+APPEND_TOKEN        
        current_val_string = current_val_string[max_char_index:]        
    hdr[key] = current_val_string

def hdr2map(hdr):
    hdr_keys = list(hdr.keys())
    my_data_map = {}
    curr_value = ''
    curr_key = ''
    for key in hdr_keys:
        separator_index = key.find(' ')
        if separator_index > 0:            
            section = key[0:separator_index]
            if not section in my_data_map:
                my_data_map[section] = {}
            ext_indx = key.find('+')
            curr_key += key[separator_index+1:ext_indx]
            val_str = str(hdr[key])
            val_last_index = val_str.find(APPEND_TOKEN)            
            curr_value += val_str[:val_last_index]            
            if ext_indx==-1:
                my_data_map[section].update({curr_key:curr_value})
                curr_value = ''
                curr_key = ''                   
    return my_data_map


class TiptopSimulation(object):
    
    def __init__(self, path, parametersFile, outputDir, outputFile, doConvolve=False,
                          doPlot=False, addSrAndFwhm=False,
                          verbose=False, getHoErrorBreakDown=False,
                          savePSDs=False):

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
        else:
            raise FileNotFoundError('No .yml or .ini can be found in '+ self.path)

        # store the parameters in data members
        self.tel_radius = self.my_data_map['telescope']['TelescopeDiameter']/2  # mas
        wvl_temp = self.my_data_map['sources_science']['Wavelength']
        if isinstance(wvl_temp, list):
            self.wvl = wvl_temp[0]  # lambda
        else:
            self.wvl = wvl_temp     # lambda
        self.zenithSrc  = self.my_data_map['sources_science']['Zenith']
        self.azimuthSrc = self.my_data_map['sources_science']['Azimuth']
        self.pointings = polarToCartesian(np.array( [self.zenithSrc, self.azimuthSrc]))
        self.xxPointigs         = self.pointings[0,:]
        self.yyPointigs         = self.pointings[1,:]
        # it checks if LO parameters are set and then it acts accordingly
        if 'sensor_LO' in self.my_data_map.keys():
            self.LOisOn = True
            if self.verbose: print('LO part is present')
        else:
            self.LOisOn = False
            self.nNaturalGS = 0
            if self.verbose: print('LO part is not present')
        if self.LOisOn:
            LO_wvl_temp = self.my_data_map['sources_LO']['Wavelength']
            if isinstance(LO_wvl_temp, list):
                self.LO_wvl = LO_wvl_temp[0]  # lambda
            else:
                self.LO_wvl = LO_wvl_temp     # lambda
            self.LO_zen     = self.my_data_map['sources_LO']['Zenith']
            self.LO_az      = self.my_data_map['sources_LO']['Azimuth']
            self.LO_fluxes  = self.my_data_map['sensor_LO']['NumberPhotons']
            self.fr         = self.my_data_map['RTC']['SensorFrameRate_LO']
            # NGSs positions
            self.NGS_flux = []
            polarNGSCoordsList = []
            for aFlux, aZen, aAz in zip(self.LO_fluxes, self.LO_zen, self.LO_az):
                polarNGSCoordsList.append([aZen, aAz])
                self.NGS_flux.append(aFlux*self.fr)
            self.polarNGSCoords     = np.asarray(polarNGSCoordsList)
            self.nNaturalGS         = self.polarNGSCoords.shape[0]

            self.cartPointingCoords = np.dstack( (self.xxPointigs, self.yyPointigs) ).reshape(-1, 2)
            cartNGSCoordsList = []
            for i in range(self.nNaturalGS):
                cartNGSCoordsList.append(polarToCartesian(self.polarNGSCoords[i,:]))
            self.cartNGSCoords = np.asarray(cartNGSCoordsList)        
        # initialize self.jitter_FWHM variable with a default value
        self.jitter_FWHM = None
        if 'jitter_FWHM' in self.my_data_map['telescope'].keys():
            self.jitter_FWHM = self.my_data_map['telescope']['jitter_FWHM']

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
        if self.addSrAndFwhm:
            for i in range(self.cubeResultsArray.shape[0]):
                hdr1['SR'+str(i).zfill(4)]   = self.getStrehl(self.cubeResultsArray[i,:,:], self.fao.ao.tel.pupil, self.fao.freq.sampRef, method='max')
            for i in range(self.cubeResultsArray.shape[0]):
                hdr1['FWHM'+str(i).zfill(4)] = getFWHM(self.cubeResultsArray[i,:,:], self.psInMas[0], method='contour', nargout=1)
            for i in range(self.cubeResultsArray.shape[0]):
                ee,rr = getEncircledEnergy(self.cubeResultsArray[i,:,:], pixelscale=self.psInMas[0], center=(self.fao.ao.cam.fovInPix/2,self.fao.ao.cam.fovInPix/2), nargout=2)
                ee_at_radius_fn = interp1d(rr, ee, kind='cubic', bounds_error=False)
                hdr1['EE50'+str(i).zfill(4)] = ee_at_radius_fn(50.0).take(0)

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
        
    def psdSetToPsfSet(self, N, freq_range, dk, mask, inputPSDs, wavelength, pixelscale, npixel, scaleFactor=1):
        NGS_SR = []
        psdArray = []
        psfLongExpArr = []
        NGS_FWHM_mas = []
        for computedPSD in inputPSDs:
            # Get the PSD at the NGSs positions at the sensing wavelength
            # computed PSD from fao are given in nm^2, i.e they are multiplied by dk**2 already
            psd            = Field(wavelength, N, freq_range, 'rad')
            psd.sampling   = computedPSD / dk**2 # the PSD must be provided in m^2.m^2
            psdArray.append(psd)
            # Get the PSF
            psfLE          = longExposurePsf(mask, psd )
            # It cuts the PSF if the PSF is larger than the requested dimension
            if psfLE.sampling.shape[0] > npixel:
                psfLE.sampling = psfLE.sampling[int(psfLE.sampling.shape[0]/2-npixel/2):int(psfLE.sampling.shape[0]/2+npixel/2),
                                                int(psfLE.sampling.shape[1]/2-npixel/2):int(psfLE.sampling.shape[1]/2+npixel/2)]
            psfLongExpArr.append(psfLE)
            # Get SR and FWHM in mas at the NGSs positions at the sensing wavelength
            s1=cpuArray(computedPSD).sum()
            SR             = np.exp(-s1*scaleFactor) # Strehl-ratio at the sensing wavelength

            NGS_SR.append(SR)
            FWHMx,FWHMy    = getFWHM( psfLE.sampling, pixelscale, method='contour', nargout=2)
            FWHM           = np.sqrt(FWHMx*FWHMy) #max(FWHMx, FWHMy) #0.5*(FWHMx+FWHMy) #average over major and minor axes

            # note : the uncertainities on the FWHM seems to create a bug in mavisLO
            NGS_FWHM_mas.append(FWHM)
            if self.verbose:
                print('SR(@',int(wavelength*1e9),'nm)  :', SR)
                print('FWHM(@',int(wavelength*1e9),'nm):', FWHM)

        return NGS_SR, psdArray, psfLongExpArr, NGS_FWHM_mas

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
        for ellp, psfLongExp in zip(self.cov_ellipses, self.psfLongExpPointingsArr):
            resSpec = residualToSpectrum(ellp, self.wvl, self.N, 1/(self.fao.ao.cam.fovInPix * self.psInMas[0]))
            self.results.append(convolve(psfLongExp, resSpec))

    def computeMetrics(self, eeRadiusInMas=50):
        if self.verbose:
            print('EE is computed for a radius of ', eeRadiusInMas,' mas')            
        self.sr, self.fwhm, self.ee = [], [], []
        for img in self.results:
            self.sr.append(getStrehl(img.sampling, self.fao.ao.tel.pupil, self.fao.freq.sampRef, method='max'))
            self.fwhm.append(getFWHM(img.sampling, self.psInMas[0], method='contour', nargout=1))
            ee_,rr_ = getEncircledEnergy(img.sampling, pixelscale=self.psInMas[0], center=(self.fao.ao.cam.fovInPix/2,self.fao.ao.cam.fovInPix/2), nargout=2)
            ee_at_radius_fn = interp1d(rr_, ee_, kind='cubic', bounds_error=False)
            self.ee.append(ee_at_radius_fn(eeRadiusInMas))

    def doOverallSimulation(self):
        self.results = []        
        # HO Part with P3 PSDs
        if self.verbose:
            print('******** HO PSD science and NGSs directions')
        self.fao = fourierModel( self.fullPathFilename, calcPSF=False, verbose=self.verbose
                           , display=False, getPSDatNGSpositions=True
                           , computeFocalAnisoCov=False, TiltFilter=self.LOisOn
                           , getErrorBreakDown=self.getHoErrorBreakDown)
        # High-order PSD caculations at the science directions and NGSs directions
        self.PSD           = self.fao.PSD # in nm^2
        self.nPointings    = self.pointings.shape[1]
        self.PSD           = self.PSD.transpose()
        self.N             = self.PSD[0].shape[0]
        self.nPixPSF       = int(self.fao.freq.nOtf /self.fao.freq.kRef_)
        self.freq_range    = self.fao.ao.cam.fovInPix*self.fao.freq.PSDstep
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
        # HO PSF
        if self.verbose:
            print('******** HO PSF')
        pointings_SR, psdPointingsArray, psfLongExpPointingsArr, pointings_FWHM_mas = self.psdSetToPsfSet(self.N, 
                                                                                                     self.freq_range, 
                                                                                                     self.dk,
                                                                                                     self.mask, 
                                                                                                     arrayP3toMastsel(self.PSD[0:self.nPointings]),
                                                                                                     self.wvl,
                                                                                                     self.psInMas[0],
                                                                                                     self.nPixPSF,
                                                                                                     scaleFactor=(2*np.pi*1e-9/self.wvl)**2)
        self.psfLongExpPointingsArr = psfLongExpPointingsArr
        # old version: if not self.LOisOn or (not self.doConvolve and not self.returnRes):
        if not self.LOisOn:
            for psfLongExp in self.psfLongExpPointingsArr:
                if self.jitter_FWHM is not None:
                    ellp = [0, sigma_from_FWHM(self.jitter_FWHM), sigma_from_FWHM(self.jitter_FWHM)]
                    self.results.append(convolve(psfLongExp,
                                   residualToSpectrum(ellp, self.wvl, N, 1/(self.fao.ao.cam.fovInPix * self.psInMas[0]))))
                else:
                    self.results.append(psfLongExp)
        else:
            # LOW ORDER PART
            if self.verbose:
                print('******** LO PART')
            self.psInMas_NGS        = self.psInMas[0] * (self.LO_wvl/self.wvl) # airy pattern PSF FWHM
            if self.verbose:
                print('******** HO PSF - NGS directions')
            NGS_SR, psdArray, psfLE_NGS, NGS_FWHM_mas = self.psdSetToPsfSet(self.N, 
                                                                       self.freq_range, 
                                                                       self.dk, 
                                                                       self.mask, 
                                                                       arrayP3toMastsel(self.PSD[-self.nNaturalGS:]), 
                                                                       self.LO_wvl, 
                                                                       self.psInMas_NGS, 
                                                                       self.nPixPSF,
                                                                       scaleFactor=(2*np.pi*1e-9/self.LO_wvl)**2)

            self.mLO           = MavisLO(self.path, self.parametersFile, verbose=self.verbose)
            self.Ctot          = self.mLO.computeTotalResidualMatrix(np.array(self.cartPointingCoords),
                                                                self.cartNGSCoords, self.NGS_flux,
                                                                NGS_SR, NGS_FWHM_mas)
            if self.doConvolve:
                self.finalConvolution()
        if self.doPlot:
            if self.LOisOn and self.doConvolve:
                tiledDisplay(self.results)
                plotEllipses(self.cartPointingCoords, self.cov_ellipses, 0.4)
            else:
                self.results[0].standardPlot(True)

        self.HO_res = np.sqrt(np.sum(self.PSD[:-self.nNaturalGS],axis=(1,2)))
        if self.LOisOn:
            self.LO_res = np.sqrt(np.trace(self.Ctot,axis1=1,axis2=2))

        self.computeOL_PSD()
        self.computeDL_PSD()
        self.cubeResults = []
        for img in self.results:
            self.cubeResults.append(img.sampling)
        self.cubeResultsArray = np.array(self.cubeResults)
                
def overallSimulation(path, parametersFile, outputDir, outputFile, doConvolve=False,
                      doPlot=False, returnRes=False, returnMetrics=False, addSrAndFwhm=False,
                      verbose=False, getHoErrorBreakDown=False, eeRadiusInMas=50,
                      savePSDs=False):
    """
        function to run the entire tiptop simulation based on the input file

        :param path2param: required, path to the parameter file.
        :type path2param: str
        :param paramFileName: required, name of the parameter file to be used without the extention.
        :type paramFileName: str
        :param outpuDir: required, path to the folder in which to write the output.
        :type outputDir: str
        :param doConvolve: optional default: False, if you want to use the natural convolution operation set to True.
        :type doConvolve: bool
        :param doPlot: optional default: False, if you want to see the result in python set this to True.
        :type doPlot: bool
        :param verbose: optional default: False, If you want all messages set this to True
        :type verbose: bool
        :param returnRes: optional default: False, The function will return the result in the environment if set to True, else it saves the result only in a .fits file.
        :type returnRes: bool
        :param returnMetrics: optional default: False, The function will return Strehl Ratio, fwhm and encircled energy within eeRadiusInMas if set to True
        :type returnMetrics: bool
        :param addSrAndFwhm: optional default: False, The function will add in the header of the fits file SR anf FWHM for each PSF.
        :type addSrAndFwhm: bool
        :param verbose: optional default: False, If you want all messages set this to True.
        :type verbose: bool
        :param getHoErrorBreakDown: optional default: False, If you want HO error breakdown set this to True.
        :type getHoErrorBreakDown: bool
        :param eeRadiusInMas: optional default: 50, used together with returnMetrics, radius used for the computation of the encirlced energy
        :type eeRadiusInMas: float
        :param savePSDs: optional default: False, If you want to save PSD in the output fits file set this to True.
        :type savePSDs: bool

        :return: TBD
        :rtype: TBD

    """
    simulation = TiptopSimulation(path, parametersFile, outputDir, outputFile, doConvolve,
                      doPlot, addSrAndFwhm, verbose, getHoErrorBreakDown, savePSDs)
    
    simulation.doOverallSimulation()

    if returnRes:
        if simulation.LOisOn:
            return simulation.HO_res, simulation.LO_res
        else:
            return simulation.HO_res  
    elif returnMetrics:
        simulation.computeMetrics(eeRadiusInMas)
        return simulation.sr, simulation.fwhm, simulation.ee
    else:
        simulation.saveResults()
