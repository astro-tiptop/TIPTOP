import os
from matplotlib import rc
from P3.aoSystem.fourierModel import *
from P3.aoSystem.FourierUtils import *
from configparser import ConfigParser

from mastsel import *

from datetime import datetime



rc("text", usetex=False)


def overallSimulation(path, parametersFile, outputDir, outputFile, doConvolve=False, doPlot=False, verbose=False, returnRes=False):
    
    # initiate the parser
    fullPathFilename = os.path.join(path, parametersFile + '.ini')    
    parser           = ConfigParser()
    parser.read(fullPathFilename);
    
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
    nPointings         = pp.shape[1]

    if verbose:
        print('******** HO PSD science and NGSs directions')
    # High-order PSD caculations at the science directions and NGSs directions
    fao = fourierModel(fullPathFilename, calcPSF=False, verbose=False, display=False, getPSDatNGSpositions=LOisOn,getErrorBreakDown=verbose)
    PSD           = fao.PSD # in nm^2
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
            #SR              = getStrehl(cp.asnumpy(psfLE.sampling), fao.ao.tel.pupil, fao.samp)
            NGS_SR.append(SR)
            FWHMx,FWHMy    = getFWHM( cp.asnumpy(psfLE.sampling),pixelscale, method='contour', nargout=2)
            FWHM           = np.sqrt(FWHMx*FWHMy) #average over major and minor axes
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
        Ctot          = mLO.computeTotalResidualMatrix(np.array(cartPointingCoords), cartNGSCoords, NGS_flux, NGS_SR, NGS_FWHM_mas)
        
        if returnRes == False:
            cov_ellipses       = mLO.ellipsesFromCovMats(Ctot)
            if verbose:
                for n in range(cov_ellipses.shape[0]): print('cov_ellipses #',n,': ',cov_ellipses[n,:], ' (unit: rad, mas, mas)')
            # FINAl CONVOLUTION
            if verbose:
                print('******** FINAl CONVOLUTION')
            results = []
            for ellp, psfLongExp in zip(cov_ellipses, psfLongExpPointingsArr):
                results.append(convolve(psfLongExp, residualToSpectrum(ellp, wvl, N, 1/(fao.ao.cam.fovInPix * fao.freq.psInMas[0]))))
    
    if doPlot:
        if LOisOn and doConvolve and returnRes == False:
            tiledDisplay(results)
            plotEllipses(cartPointingCoords, cov_ellipses, 0.4)
        else:
            tiledDisplay(results)
    
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