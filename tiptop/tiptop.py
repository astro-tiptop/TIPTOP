import os
import numpy as np
from matplotlib import rc
from P3.aoSystem.fourierModel import *
from P3.aoSystem.FourierUtils import *
from configparser import ConfigParser
import yaml

from mastsel import *

from datetime import datetime

rc("text", usetex=False)

def overallSimulation(path, parametersFile, outputDir, outputFile, doConvolve=False, 
                      doPlot=False, verbose=False, returnRes=False):
    """
    function to run the entire tiptop simulation based on the imput file
    
    :param path2param: required, path to the parameter file
    :type path2param: str
    :param paramFileName: required, name of the parameter file to be used without the extention
    :type paramFileName: str
    :param outpuDir: required, path to the folder in which to write the output
    :type outputDir: str
    :param doConvolve: optional default: False, if you want to use the natural convolution operation  set to True
    :type doConvolve: bool
    :param doPlot: optional default: False, if you want to see the result in python set this to True
    :type doPlot: bool
    :param verbose: optional default: False, If you want all messages set this to True
    :type verbose: bool
    :param returnRes: optionnal default: False, The function will return the result in the environment if set to True, else it saves the result only in a .fits file. 

    :return: TBD
    :rtype: TBD

    """
    #TODO remove this prints in one of the next releases of the library
    print('ATTENTION: interface of this function is changed.')
    print('           windPsdFile is no more an input of overallSimulation,')
    print('           it must be set as a parameter in the telescope section of the ini file.')
          

    # initiate the parser
    fullPathFilename_ini = os.path.join(path, parametersFile + '.ini')
    fullPathFilename_yml = os.path.join(path, parametersFile + '.yml')

    if os.path.exists(fullPathFilename_yml):
        with open(fullPathFilename_yml) as f:
            my_yaml_dict = yaml.safe_load(f)
        # read main parameters
        tel_radius = my_yaml_dict['telescope']['TelescopeDiameter']/2  # mas
        wvl_temp = my_yaml_dict['sources_science']['Wavelength']
        if isinstance(wvl_temp, list):
            wvl = wvl_temp[0]  # lambda
        else:
            wvl = wvl_temp     # lambda
        zenithSrc  = my_yaml_dict['sources_science']['Zenith']
        azimuthSrc = my_yaml_dict['sources_science']['Azimuth']
        
        # it checks if LO parameters are set and then it acts accordingly
        if 'sensor_LO' in my_yaml_dict.keys():
            LOisOn = True
            if verbose: print('LO part is present')
        else:
            LOisOn = False
            nNaturalGS = 0
            if verbose: print('LO part is not present')

        if LOisOn: 
            LO_wvl_temp = my_yaml_dict['sources_LO']['Wavelength']
            if isinstance(LO_wvl_temp, list):
                LO_wvl = LO_wvl_temp[0]  # lambda
            else:
                LO_wvl = LO_wvl_temp     # lambda
            LO_zen     = my_yaml_dict['sources_LO']['Zenith']
            LO_az      = my_yaml_dict['sources_LO']['Azimuth']
            LO_fluxes  = my_yaml_dict['sensor_LO']['NumberPhotons']
            fr         = my_yaml_dict['RTC']['SensorFrameRate_LO']

        fao = fourierModel( fullPathFilename_yml, calcPSF=False, verbose=False
                           , display=False, getPSDatNGSpositions=True)
        
    elif os.path.exists(fullPathFilename_ini):
        parser           = ConfigParser()
        parser.read(fullPathFilename_ini);
        # read main parameters
        tel_radius = eval(parser.get('telescope', 'TelescopeDiameter'))/2  # mas
        wvl_temp = eval(parser.get('sources_science', 'Wavelength'))
        if isinstance(wvl_temp, list):
            wvl = wvl_temp[0]  # lambda
        else:
            wvl = wvl_temp     # lambda
        zenithSrc  = eval(parser.get('sources_science', 'Zenith'))
        azimuthSrc = eval(parser.get('sources_science', 'Azimuth'))
        # it checks if LO parameters are set and then it acts accordingly
        if parser.has_section('sensor_LO'):
            LOisOn = True
            if verbose: print('LO part is present')
        else:
            LOisOn = False
            nNaturalGS = 0
            if verbose: print('LO part is not present')

        if LOisOn:
            LO_wvl_temp = eval(parser.get('sources_LO', 'Wavelength'))
            if isinstance(LO_wvl_temp, list):
                LO_wvl = LO_wvl_temp[0]  # lambda
            else:
                LO_wvl = LO_wvl_temp     # lambda
            LO_zen     = eval(parser.get('sources_LO', 'Zenith')) 
            LO_az      = eval(parser.get('sources_LO', 'Azimuth'))
            LO_fluxes  = eval(parser.get('sensor_LO', 'NumberPhotons'))
            fr         = eval(parser.get('RTC', 'SensorFrameRate_LO'))


        fao = fourierModel( fullPathFilename_ini, calcPSF=False, verbose=False
                           , display=False, getPSDatNGSpositions=True)
    else:
        raise FileNotFoundError('No .yml or .ini can be found in '+ path)

    if LOisOn:        
        # NGSs positions
        NGS_flux = []
        polarNGSCoordsList = []
        for aFlux, aZen, aAz in zip(LO_fluxes, LO_zen, LO_az):
            polarNGSCoordsList.append([aZen, aAz])   
            NGS_flux.append(LO_fluxes[0]*fr)
            polarNGSCoords     = np.asarray(polarNGSCoordsList)
            nNaturalGS         = polarNGSCoords.shape[0]

    pp                 = polarToCartesian(np.array( [zenithSrc, azimuthSrc]))
    xxPointigs         = pp[0,:]
    yyPointigs         = pp[1,:]
        
#  cartPointingCoords=np.array([xxPointigs, yyPointigs]).transpose(),
#  extraPSFsDirections=polarNGSCoordsList
#  kcExt=kcExt
#  pitchScaling=pitchScaling
#  path_pupil=path_pupil

    # High-order PSD caculations at the science directions and NGSs directions
#    PSD           = fao.powerSpectrumDensity() # in nm^2 : old way?
    PSD           = fao.PSD # in nm^2
    nPointings         = pp.shape[1]
    if verbose:
        print('******** HO PSD science and NGSs directions')

    PSD           = PSD.transpose()
    N             = PSD[0].shape[0]
    freq_range    = fao.ao.cam.fovInPix*fao.freq.PSDstep 
    pitch         = 1/freq_range
    grid_diameter = pitch*N
    sx            = int(2*np.round(tel_radius/pitch))
    dk            = 1e9*fao.freq.kcMax_/fao.freq.resAO
    
    # Define the pupil shape
    mask = Field(wvl, N, grid_diameter)
    mask.sampling = cp.asarray(congrid(fao.ao.tel.pupil, [sx, sx]))
    mask.sampling = zeroPad(mask.sampling, (N-sx)//2)
    if verbose:
        print('fao.samp:', fao.freq.samp)
        print('fao.freq.psInMas:', fao.freq.psInMas)

    def psdSetToPsfSet(inputPSDs,wavelength,pixelscale,scaleFactor=1,verbose=False):
        NGS_SR = []
        psdArray = []
        psfLongExpArr = []
        NGS_FWHM_mas = []
        for computedPSD in inputPSDs:
            # Get the PSD at the NGSs positions at the sensing wavelength
            # computed PSD from fao are given in nm^2, i.e they are multiplied by dk**2 already
            psd            = Field(wavelength, N, freq_range, 'rad')
            psd.sampling   = cp.asarray( computedPSD / dk**2) # the PSD must be provided in m^2.m^2
            psdArray.append(psd)
            # Get the PSF
            psfLE          = longExposurePsf(mask, psd )     
            psfLongExpArr.append(psfLE)
            # Get SR and FWHM in mas at the NGSs positions at the sensing wavelength
            SR             = np.exp(-computedPSD.sum()* scaleFactor) # Strehl-ratio at the sensing wavelength
            NGS_SR.append(SR)
            FWHMx,FWHMy    = getFWHM( psfLE.sampling, pixelscale, method='contour', nargout=2)
            FWHM           = max(FWHMx, FWHMy) #0.5*(FWHMx+FWHMy) #average over major and minor axes
            # note : the uncertainities on the FWHM seems to create a bug in mavisLO
            NGS_FWHM_mas.append(FWHM)
            if verbose:
                print('SR(@',int(wavelength*1e9),'nm)  :', SR)
                print('FWHM(@',int(wavelength*1e9),'nm):', FWHM)
            
        return NGS_SR, psdArray, psfLongExpArr, NGS_FWHM_mas

    # HO PSF
    if verbose:
        print('******** HO PSF')
    pointings_SR, psdPointingsArray, psfLongExpPointingsArr, pointings_FWHM_mas = psdSetToPsfSet(PSD[0:nPointings],
                                                                                                 wvl,
                                                                                                 fao.freq.psInMas[0],
                                                                                                 scaleFactor=(2*np.pi*1e-9/wvl)**2,
                                                                                                 verbose=verbose)
    
    if LOisOn == False or (doConvolve == False and returnRes == False):
        results = []
        for psfLongExp in psfLongExpPointingsArr:
            results.append(psfLongExp)
    else:
        # LOW ORDER PART
        if verbose:
            print('******** LO PART')
        psInMas_NGS        = fao.freq.psInMas[0] * (LO_wvl/wvl) #airy pattern PSF FWHM
        
        if verbose:
            print('******** HO PSF - NGS directions')
        NGS_SR, psdArray, psfLE_NGS, NGS_FWHM_mas = psdSetToPsfSet(PSD[-nNaturalGS:],
                                                                   LO_wvl, psInMas_NGS,
                                                                   scaleFactor=(2*np.pi*1e-9/LO_wvl)**2,
                                                                   verbose=verbose)
        
        cartPointingCoords = np.dstack( (xxPointigs, yyPointigs) ).reshape(-1, 2)
        cartNGSCoordsList = []
        for i in range(nNaturalGS):
            cartNGSCoordsList.append(polarToCartesian(polarNGSCoords[i,:]))
        cartNGSCoords = np.asarray(cartNGSCoordsList)
        mLO           = MavisLO(path, parametersFile,verbose=verbose)
        Ctot          = mLO.computeTotalResidualMatrix(np.array(cartPointingCoords)
                                                       , cartNGSCoords, NGS_flux, 
                                                       NGS_SR, NGS_FWHM_mas)
        
        if returnRes == False:
            cov_ellipses       = mLO.ellipsesFromCovMats(Ctot)
            if verbose:
                for n in range(cov_ellipses.shape[0]): 
                    print('cov_ellipses #',n,': ',cov_ellipses[n,:], ' (unit: rad, mas, mas)')
            # FINAl CONVOLUTION
            if verbose:
                print('******** FINAl CONVOLUTION')
            results = []
            for ellp, psfLongExp in zip(cov_ellipses, psfLongExpPointingsArr):
                results.append(convolve(psfLongExp, residualToSpectrum(ellp, wvl
                                                                       , N, 1/(fao.ao.cam.fovInPix 
                                                                               * fao.freq.psInMas[0]))))
    
    if doPlot:
        if LOisOn and doConvolve and returnRes == False:
            tiledDisplay(results)
            plotEllipses(cartPointingCoords, cov_ellipses, 0.4)
        else:
            results[0].standardPlot(True)

    # save PSF cube in fits    
    hdul1 = fits.HDUList()
    cube =[]
    hdul1.append(fits.PrimaryHDU())
    for img in results:
        cube.append(img.sampling)
    
    if returnRes:
        HO_res = np.sqrt(np.sum(PSD[:-nNaturalGS],axis=(1,2)))
        if LOisOn:
            LO_res = np.sqrt(np.trace(Ctot,axis1=1,axis2=2))
    
            return HO_res, LO_res
        else:
            return HO_res
    else:
        # save PSF cube in fits
        hdul1 = fits.HDUList()
        cube =[]
        hdul1.append(fits.PrimaryHDU())
        for img in results:
            cube.append(cp.asnumpy(img.sampling))
        
        cube = np.array(cube)
        hdul1.append(fits.ImageHDU(data=cube))
        hdul1.append(fits.ImageHDU(data=pp)) # append cartesian coordinates
    
        #############################
        # header
        hdr0 = hdul1[0].header
        now = datetime.now()
        hdr0['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        # header of the PSFs
        hdr1 = hdul1[1].header
        hdr1['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        hdr1['CONTENT'] = "PSF CUBE"
        hdr1['SIZE'] = str(cube.shape)
        hdr1['WL_NM'] = str(int(wvl*1e9))
        hdr1['PIX_MAS'] = str(fao.freq.psInMas[0])
        # header of the coordinates
        hdr2 = hdul1[2].header
        hdr2['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        hdr2['CONTENT'] = "CARTESIAN COORD. IN ASEC OF THE SOURCES"
        hdr2['SIZE'] = str(pp.shape)
        #############################

        hdul1.writeto( os.path.join(outputDir, outputFile + '.fits'), overwrite=True)
        print("Output cube shape:", cube.shape)
