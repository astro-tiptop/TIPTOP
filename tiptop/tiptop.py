import os
import numpy as np
from matplotlib import rc
from p3.aoSystem.fourierModel import *
from p3.aoSystem.FourierUtils import *
from configparser import ConfigParser
import yaml

from mastsel import *

from datetime import datetime

rc("text", usetex=False)

def overallSimulation(path, parametersFile, outputDir, outputFile, doConvolve=False,
                      doPlot=False, returnRes=False, addSrAndFwhm=False,
                      verbose=False, getHoErrorBreakDown=False):
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
    :param returnRes: optional default: False, The function will return the result in the environment if set to True, else it saves the result only in a .fits file.
    :type returnRes: bool
    :param addSrAndFwhm: optional default: False, The function will add in the header of the fits file SR anf FWHM for each PSF.
    :type addSrAndFwhm: bool
    :param verbose: optional default: False, If you want all messages set this to True
    :type verbose: bool
    :param getHoErrorBreakDown: optional default: False, If you want HO error breakdosn set this to True
    :type getHoErrorBreakDown: bool

    :return: TBD
    :rtype: TBD

    """

    # initiate the parser
    fullPathFilename_ini = os.path.join(path, parametersFile + '.ini')
    fullPathFilename_yml = os.path.join(path, parametersFile + '.yml')

    # initialize jitter_FWHM variable with a default value
    jitter_FWHM = None

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

        if 'jitter_FWHM' in my_yaml_dict['telescope'].keys():
            jitter_FWHM = my_yaml_dict['telescope']['jitter_FWHM']

        fao = fourierModel( fullPathFilename_yml, calcPSF=False, verbose=verbose
                           , display=False, getPSDatNGSpositions=True
                           , computeFocalAnisoCov=False, TiltFilter=LOisOn
                           , getErrorBreakDown=getHoErrorBreakDown)

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
        if parser.has_option('telescope', 'jitter_FWHM'):
            jitter_FWHM = eval(parser.get('telescope', 'jitter_FWHM'))

        fao = fourierModel( fullPathFilename_ini, calcPSF=False, verbose=verbose
                           , display=False, getPSDatNGSpositions=True
                           , computeFocalAnisoCov=False, TiltFilter=LOisOn
                           , getErrorBreakDown=getHoErrorBreakDown)
        
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

    # High-order PSD caculations at the science directions and NGSs directions
    PSD                = fao.PSD # in nm^2
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
            FWHM           = np.sqrt(FWHMx*FWHMy) #max(FWHMx, FWHMy) #0.5*(FWHMx+FWHMy) #average over major and minor axes
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
            if jitter_FWHM is not None:
                ellp = [0, sigma_from_FWHM(jitter_FWHM), sigma_from_FWHM(jitter_FWHM)]
                results.append(convolve(psfLongExp,
                               residualToSpectrum(ellp, wvl, N, 1/(fao.ao.cam.fovInPix * fao.freq.psInMas[0]))))
            else:
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
                print('******** FINAL CONVOLUTION')
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
        # OPEN-LOOP PSD
        k   = np.sqrt(fao.freq.k2_)
        pf  = FourierUtils.pistonFilter(fao.ao.tel.D,k)
        psdOL = Field(wvl, N, freq_range, 'rad')
        temp = fao.ao.atm.spectrum(k) * pf
        psdOL.sampling = cp.asarray(fao.ao.atm.spectrum(k) * pf * (fao.freq.wvlRef/np.pi)**2) # the PSD must be provided in m^2.m^2
        # Get the OPEN-LOOP PSF
        psfOL = longExposurePsf(mask, psdOL)
        if doPlot:
            fig, ax1 = plt.subplots(1,1)
            im = ax1.imshow(np.log(np.abs(psfOL.sampling) + 1e-20), cmap='hot')
            ax1.set_title('open loop PSF', color='black')
            
        # DIFFRACTION LIMITED PSD
        psdDL = Field(wvl, N, freq_range, 'rad')
        psfDL = longExposurePsf(mask, psdDL)
        if doPlot:
            fig, ax2 = plt.subplots(1,1)
            im = ax2.imshow(np.log(np.abs(psfDL.sampling) + 1e-20), cmap='hot')
            ax2.set_title('diffraction limited PSF', color='black')

        # save PSF cube in fits
        hdul1 = fits.HDUList()
        cube =[]
        hdul1.append(fits.PrimaryHDU())
        for img in results:
            cube.append(img.sampling)

        cube = np.array(cube)
        hdul1.append(fits.ImageHDU(data=cube))
        hdul1.append(fits.ImageHDU(data=psfOL.sampling)) # append open-loop PSF
        hdul1.append(fits.ImageHDU(data=psfDL.sampling)) # append diffraction limited PSF

        #############################
        # header
        hdr0 = hdul1[0].header
        now = datetime.now()
        hdr0['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        # parameters in the header
        if os.path.exists(fullPathFilename_yml):
            for key_primary in my_yaml_dict:
                for key_secondary in my_yaml_dict[key_primary]:
                    temp = my_yaml_dict[key_primary][key_secondary]
                    if isinstance(temp, list):
                        iii = 0
                        for elem in temp:
                            hdr0['HIERARCH '+key_primary+' '+key_secondary +' '+str(iii)] = elem
                            iii += 1
                    else:
                        hdr0['HIERARCH '+key_primary+' '+key_secondary] = temp
        # header of the PSFs
        hdr1 = hdul1[1].header
        hdr1['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        hdr1['CONTENT'] = "PSF CUBE"
        hdr1['SIZE'] = str(cube.shape)
        hdr1['WL_NM'] = str(int(wvl*1e9))
        hdr1['PIX_MAS'] = str(fao.freq.psInMas[0])
        hdr1['CC'] = "CARTESIAN COORD. IN ASEC OF THE "+str(pp.shape[1])+" SOURCES"
        for i in range(pp.shape[1]):
            hdr1['CCX'+str(i).zfill(4)] = pp[0,i]
            hdr1['CCY'+str(i).zfill(4)] = pp[1,i]
        if addSrAndFwhm:
            for i in range(cube.shape[0]):
                hdr1['SR'+str(i).zfill(4)]   = getStrehl(cube[i,:,:], fao.ao.tel.pupil, fao.freq.sampRef)
            for i in range(cube.shape[0]):
                hdr1['FWHM'+str(i).zfill(4)] = getFWHM(cube[i,:,:], fao.freq.psInMas[0], method='contour', nargout=1)

        # header of the OPEN-LOOP PSF
        hdr2 = hdul1[2].header
        hdr2['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        hdr2['CONTENT'] = "OPEN-LOOP PSF"
        hdr2['SIZE'] = str(psfOL.sampling.shape)
        
        # header of the DIFFRACTION LIMITED PSF
        hdr3 = hdul1[3].header
        hdr3['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        hdr3['CONTENT'] = "DIFFRACTION LIMITED PSF"
        hdr3['SIZE'] = str(psfDL.sampling.shape)
        #############################

        hdul1.writeto( os.path.join(outputDir, outputFile + '.fits'), overwrite=True)
        if verbose:
            print("Output cube shape:", cube.shape)
