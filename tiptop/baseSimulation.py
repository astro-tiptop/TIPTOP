from p3.aoSystem.fourierModel import *
from p3.aoSystem.FourierUtils import *
from mastsel import *

from .tiptopUtils import *

from matplotlib import cm
import matplotlib as mpl
norm = mpl.colors.Normalize(vmin=0, vmax=1)
rc("text", usetex=False)


class baseSimulation(object):
    
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


    def configLO(self, astIndex=None):
        self.cartSciencePointingCoords = np.dstack( (self.xxSciencePointigs, self.yySciencePointigs) ).reshape(-1, 2)
        self.fr         = self.my_data_map['RTC']['SensorFrameRate_LO']
        # Here we assume the same wavelenght for all the phon counts of the stars in the asterism
        LO_wvl_temp = self.my_data_map['sources_LO']['Wavelength']

        if isinstance(LO_wvl_temp, list):
            self.LO_wvl = LO_wvl_temp[0]  # lambda
        else:
            self.LO_wvl = LO_wvl_temp     # lambda

        self.LO_zen_field    = self.my_data_map['sources_LO']['Zenith']
        self.LO_az_field     = self.my_data_map['sources_LO']['Azimuth']
        self.LO_fluxes_field = self.my_data_map['sensor_LO']['NumberPhotons']

        self.NGS_fluxes_field = []
        polarNGSCoordsList = []
        for aFlux, aZen, aAz in zip(self.LO_fluxes_field, self.LO_zen_field, self.LO_az_field):
            polarNGSCoordsList.append([aZen, aAz])
            self.NGS_fluxes_field.append(aFlux*self.fr)
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
        self.NGS_fluxes_asterism = []
        for iid in self.currentAsterismIndices:
            self.NGS_fluxes_asterism.append(self.NGS_fluxes_field[iid])
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
        if self.addSrAndFwhm:
            for i in range(self.cubeResultsArray.shape[0]):
                hdr1['SR'+str(i).zfill(4)]   = float(getStrehl(self.cubeResultsArray[i,:,:], self.fao.ao.tel.pupil, self.fao.freq.sampRef, method='max'))
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


    def psdSetToPsfSet(self, N, freq_range, dk, mask, inputPSDs, wavelength, pixelscale, npixel, scaleFactor=1, oversampling=1):
        sources_SR = []
        psdArray = []
        psfLongExpArr = []
        sources_FWHM_mas = []

        for computedPSD in inputPSDs:
            # Get the PSD at the NGSs positions at the sensing wavelength
            # computed PSD from fao are given in nm^2, i.e they are multiplied by dk**2 already
            psd            = Field(wavelength, N, freq_range, 'rad')
            psd.sampling   = computedPSD / dk**2 # the PSD must be provided in m^2.m^2
            psdArray.append(psd)
            # Get the PSF
            psfLE          = longExposurePsf(mask, psd )
            
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
            s1=cpuArray(computedPSD).sum()
            SR             = np.exp(-s1*scaleFactor) # Strehl-ratio at the sensing wavelength

            sources_SR.append(SR)
            FWHMx,FWHMy    = getFWHM( psfLE.sampling, pixelscale, method='contour', nargout=2)
            FWHM           = np.sqrt(FWHMx*FWHMy) #max(FWHMx, FWHMy) #0.5*(FWHMx+FWHMy) #average over major and minor axes

            # note : the uncertainities on the FWHM seems to create a bug in mavisLO
            sources_FWHM_mas.append(FWHM)
            if self.verbose:
                print('SR(@',int(wavelength*1e9),'nm)  :', SR)
                print('FWHM(@',int(wavelength*1e9),'nm):', FWHM)

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


    def doOverallSimulation(self, astIndex=None):

        if self.LOisOn:        
            self.configLO(astIndex)
        
        self.results = []
        # HO Part with P3 PSDs
        if self.verbose:
            print('******** HO PSD science and NGSs directions')
        if astIndex is None or astIndex==0:
            self.fao = fourierModel( self.fullPathFilename, calcPSF=False, verbose=self.verbose
                               , display=False, getPSDatNGSpositions=True
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
                                                                                                         scaleFactor=(2*np.pi*1e-9/self.wvl)**2,
                                                                                                         oversampling=self.oversampling)
            self.psfLongExpPointingsArr = psfLongExpPointingsArr
            
        # old version: if not self.LOisOn or (not self.doConvolve and not self.returnRes):
        if not self.LOisOn:
            for psfLongExp in self.psfLongExpPointingsArr:
                if self.jitter_FWHM is not None:
                    ellp = [0, sigma_from_FWHM(self.jitter_FWHM), sigma_from_FWHM(self.jitter_FWHM)]
                    self.results.append(convolve(psfLongExp,
                                   residualToSpectrum(ellp, self.wvl, self.nPixPSF, 1/(self.fao.ao.cam.fovInPix * self.psInMas[0]))))
                else:
                    self.results.append(psfLongExp)
        else:
            if astIndex is None or astIndex==0:
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
                                                                           arrayP3toMastsel(self.PSD[-self.nNaturalGS_field:]), 
                                                                           self.LO_wvl, 
                                                                           self.psInMas_NGS, 
                                                                           self.nPixPSF,
                                                                           scaleFactor=(2*np.pi*1e-9/self.LO_wvl)**2,
                                                                           oversampling=self.oversampling)
                self.NGS_SR_field           = NGS_SR
                self.NGS_FWHM_mas_field     = NGS_FWHM_mas
                self.mLO           = MavisLO(self.path, self.parametersFile, verbose=self.verbose)

            if astIndex is None:
                self.Ctot          = self.mLO.computeTotalResidualMatrix(np.array(self.cartSciencePointingCoords),
                                                                         self.cartNGSCoords_field, self.NGS_fluxes_field,
                                                                         self.NGS_SR_field, self.NGS_FWHM_mas_field, doAll=True)
            else:
                self.NGS_SR_asterism = []
                for iid in self.currentAsterismIndices:
                    self.NGS_SR_asterism.append(self.NGS_SR_field[iid])
                self.NGS_FWHM_mas_asterism = []
                for iid in self.currentAsterismIndices:
                    self.NGS_FWHM_mas_asterism.append(self.NGS_FWHM_mas_field[iid])
                if astIndex==0:
                    self.mLO.computeTotalResidualMatrix(np.array(self.cartSciencePointingCoords),
                                                         self.cartNGSCoords_field, self.NGS_fluxes_field,
                                                         self.NGS_SR_field, self.NGS_FWHM_mas_field, doAll=False)
                self.Ctot          = self.mLO.computeTotalResidualMatrixI(self.currentAsterismIndices,
                                                                          np.array(self.cartSciencePointingCoords),
                                                                          np.array(self.cartNGSCoords_asterism), self.NGS_fluxes_asterism,
                                                                          self.NGS_SR_asterism, self.NGS_FWHM_mas_asterism)
                
            if self.doConvolve:
                if astIndex is None:
                    self.finalConvolution()
                else:
                    self.cov_ellipses = self.mLO.ellipsesFromCovMats(self.Ctot)
 
        if self.doPlot:
            if self.LOisOn and self.doConvolve:
                tiledDisplay(self.results)
                plotEllipses(self.cartSciencePointingCoords, self.cov_ellipses, 0.4)
            else:
                self.results[0].standardPlot(True)

        if astIndex is None:
            self.HO_res = np.sqrt(np.sum(self.PSD[:-self.nNaturalGS_field],axis=(1,2)))
            if self.LOisOn:
                self.LO_res = np.sqrt(np.trace(self.Ctot,axis1=1,axis2=2))
            self.computeOL_PSD()
            self.computeDL_PSD()
            self.cubeResults = []
            for img in self.results:
                self.cubeResults.append(img.sampling)
            self.cubeResultsArray = np.array(self.cubeResults)
