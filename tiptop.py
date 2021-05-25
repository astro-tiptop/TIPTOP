import os
from matplotlib import rc
from P3.aoSystem.fourierModel import *
from P3.aoSystem.FourierUtils import *
from configparser import ConfigParser

from mavis import *

rc("text", usetex=False)

def overallSimulation(path, parametersFile, windPsdFile, outputDir, outputFile, path_pupil='', pitchScaling=1,doConvolve=False, doPlot=False):
    
    # initiate the parser
    fullPathFilename = os.path.join(path, parametersFile + '.ini')    
    parser           = ConfigParser()
    parser.read(fullPathFilename);
    
    # read main parameters
    tel_radius = eval(parser.get('telescope', 'TelescopeDiameter'))/2      # mas
    wvl        = eval(parser.get('sources_science', 'Wavelength'))[0]  # lambda
    zenithSrc  = eval(parser.get('sources_science', 'Zenith'))
    azimuthSrc = eval(parser.get('sources_science', 'Azimuth'))
    wvl_LO     = eval(parser.get('sources_LO', 'Wavelength'))  # lambda
    LO_zen     = eval(parser.get('sources_LO', 'Zenith')) 
    LO_az      = eval(parser.get('sources_LO', 'Azimuth'))
    fluxes     = eval(parser.get('sensor_LO', 'NumberPhotons'))
    fr         = eval(parser.get('RTC', 'SensorFrameRate_LO'))
    
    # NGSs positions
    NGS_flux = []
    polarNGSCoordsList = []
    for aFlux, aZen, aAz in zip(fluxes, LO_zen, LO_az):
        polarNGSCoordsList.append([aZen, aAz])   
        NGS_flux.append(fluxes[0]*fr)

    polarNGSCoords     = np.asarray(polarNGSCoordsList)
    nNaturalGS         = polarNGSCoords.shape[0]

    pp                 = polarToCartesian(np.array( [zenithSrc, azimuthSrc]))
    xxPointigs         = pp[0,:]
    yyPointigs         = pp[1,:]

    # bypass the dm cut off frequency definition
    if pitchScaling != 1:
        pitchs_dm    = np.array(eval(parser.get('DM', 'DmPitchs')))
        kcExt        = 1/(2*pitchs_dm.min()*pitchScaling)
    else:
        kcExt = None
        
#  cartPointingCoords=np.array([xxPointigs, yyPointigs]).transpose(),
#  extraPSFsDirections=polarNGSCoordsList
#  kcExt=kcExt
#  pitchScaling=pitchScaling
#  path_pupil=path_pupil

    # High-order PSD caculations at the science directions and NGSs directions
    fao = fourierModel(fullPathFilename, calcPSF=False, verbose=False, display=False, getPSDatNGSpositions=True)
    PSD           = fao.powerSpectrumDensity() # in nm^2
    PSD           = PSD.transpose()
    N             = PSD[0].shape[0]
    freq_range    = fao.ao.wfs.detector.fovInPix*fao.freq.PSDstep # fao.psf_FoV/fao.wvlRef/206264.8
    pitch         = 1/freq_range
    grid_diameter = pitch*N
    sx            = int(2*np.round(tel_radius/pitch))
    dk            = 1e9*fao.freq.kc_/fao.freq.resAO
    
    # Define the pupil shape
    mask = Field(wvl, N, grid_diameter)
    mask.sampling = cp.asarray(congrid(fao.ao.tel.pupil, [sx, sx]))
    mask.sampling = zeroPad(mask.sampling, (N-sx)//2)
    #mask.standardPlot(False)
    print('fao.samp:', fao.freq.samp)
    
    def psdSetToPsfSet(inputPSDs,wavelength,pixelscale,scaleFactor=1,verbose=False):
        NGS_SR = []
        psdArray = []
        psfLongExpArr = []
        NGS_FWHM_mas = []
        for computedPSD in inputPSDs:    
            # Get the PSD at the NGSs positions at the sensing wavelength
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
            FWHM           = 0.5*(FWHMx+FWHMy) #average over major and minor axes
            NGS_FWHM_mas.append(FWHM)
            if verbose == True:
                print('SR at the sensing wavelength:', SR)            
                print('FWHM at the sensing wavelength:', FWHM)            
            
        return NGS_SR, psdArray, psfLongExpArr, NGS_FWHM_mas

    # HO PSF
    pointings_SR, psdPointingsArray, psfLongExpPointingsArr, pointings_FWHM_mas = psdSetToPsfSet(PSD[:-nNaturalGS],wvl,fao.freq.psInMas,scaleFactor=(2*np.pi*1e-9/wvl)**2, verbose=True)
    
    if doConvolve == False:
        results = psfLongExpPointingsArr
    else:
        # LOW ORDER PART
        psInMas_NGS        = fao.freq.psInMas * (wvl_LO/wvl) #airy pattern PSF FWHM
        NGS_SR, psdArray, psfLE_NGS, NGS_FWHM_mas = psdSetToPsfSet(PSD[-nNaturalGS:],wvl_LO,psInMas_NGS,scaleFactor=(2*np.pi*1e-9/wvl_LO)**2)
        cartPointingCoords = np.dstack( (xxPointigs, yyPointigs) ).reshape(-1, 2)
        cartNGSCoordsList = []
        for i in range(nNaturalGS):
            cartNGSCoordsList.append(polarToCartesian(polarNGSCoords[i,:]))        
        cartNGSCoords = np.asarray(cartNGSCoordsList)
        mLO                = MavisLO(path, parametersFile, windPsdFile)
        Ctot               = mLO.computeTotalResidualMatrix(np.array(cartPointingCoords), cartNGSCoords, NGS_flux, NGS_SR, NGS_FWHM_mas)
        cov_ellipses       = mLO.ellipsesFromCovMats(Ctot)
        # FINAl CONVOLUTION
        results = []
        for ellp, psfLongExp in zip(cov_ellipses, psfLongExpPointingsArr):
            results.append(convolve(psfLongExp, residualToSpectrum(ellp, wvl, N, 1/(fao.ao.wfs.detector.fovInPix * fao.freq.psInMas))))

    if doPlot:
        if doConvolve:
            tiledDisplay(results)
            plotEllipses(cartPointingCoords, cov_ellipses, 0.4)
        else:
            results[0].standardPlot(True)

    # save PSF cube in fits    
    hdul1 = fits.HDUList()
    cube =[]
    hdul1.append(fits.PrimaryHDU())
    for img in results:
        cube.append(cp.asnumpy(img.sampling))
    
    cube = np.array(cube)
    hdul1.append(fits.ImageHDU(data=cube))
    hdul1.writeto( os.path.join(outputDir, outputFile + '.fits'), overwrite=True)
    print("Output cube shape:", cube.shape)