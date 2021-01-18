from matplotlib import rc
from fourierPSF.fourierModel import *
from fourierPSF.FourierUtils import *
from configparser import ConfigParser

from mavis import *

rc("text", usetex=False)

def overallSimulation(path, parametersFile, windPsdFile, outputDir, outputFile, path_pupil='', pitchScaling=1,doConvolve=False):
    
    # initiate the parser
    fullPathFilename = path + parametersFile + '.ini'    
    parser           = ConfigParser()
    parser.read(fullPathFilename);
    
    # read main parameters
    wvl        = eval(parser.get('PSF_DIRECTIONS', 'ScienceWavelength'))[0]  # lambda
    wvl_LO     = eval(parser.get('SENSOR_LO', 'SensingWavelength_LO'))  # lambda
    pixel_psf  = eval(parser.get('PSF_DIRECTIONS', 'psInMas'))
    tel_radius = eval(parser.get('telescope', 'TelescopeDiameter'))/2      # mas
    LO_zen     = eval(parser.get('SENSOR_LO', 'GuideStarZenith_LO')) 
    LO_az      = eval(parser.get('SENSOR_LO', 'GuideStarAzimuth_LO'))
    fluxes     = eval(parser.get('SENSOR_LO', 'nph_LO'))
    fr         = eval(parser.get('SENSOR_LO', 'SensorFrameRate_LO'))
    zenithSrc  = np.array(eval(parser.get('PSF_DIRECTIONS', 'ScienceZenith')))
    azimuthSrc = np.array(eval(parser.get('PSF_DIRECTIONS', 'ScienceAzimuth')))
    
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
        
    
    # High-order PSD caculations at the science directions and NGSs directions
    fao = fourierModel(fullPathFilename, calcPSF=False, verbose=False, display=False, \
                   cartPointingCoords=np.array([xxPointigs, yyPointigs]).transpose() , extraPSFsDirections=polarNGSCoordsList, \
                      kcExt=kcExt,pitchScaling=pitchScaling,path_pupil=path_pupil)
    PSD           = fao.powerSpectrumDensity() # in nm^2
    PSD           = PSD.transpose()
    N             = PSD[0].shape[0]    
    freq_range    = fao.fovInPixel*fao.PSDstep # fao.psf_FoV/fao.wvlRef/206264.8
    pitch         = 1/freq_range
    grid_diameter = pitch*N
    sx            = int(2*np.round(tel_radius/pitch))
    dk            = 1e9*fao.kc/fao.resAO
    
    # Define the pupil shape
    mask = Field(wvl, N, grid_diameter)
    mask.sampling = cp.asarray(congrid(fao.tel.pupil, [sx, sx]))
    mask.sampling = zeroPad(mask.sampling, (N-sx)//2)
    #mask.standardPlot(False)
    print('fao.samp:', fao.samp)
    
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
            #SR              = getStrehl(cp.asnumpy(psfLE.sampling), fao.tel.pupil, fao.samp)
            NGS_SR.append(SR)
            FWHMx,FWHMy    = getFWHM( cp.asnumpy(psfLE.sampling),pixelscale, method='contour', nargout=2)
            FWHM           = 0.5*(FWHMx+FWHMy) #average over major and minor axes
            NGS_FWHM_mas.append(FWHM)
            if verbose == True:
                print('SR at the sensing wavelength:', SR)            
                print('FWHM at the sensing wavelength:', FWHM)            
            
        return NGS_SR, psdArray, psfLongExpArr, NGS_FWHM_mas

    # HO PSF
    pointings_SR, psdPointingsArray, psfLongExpPointingsArr, pointings_FWHM_mas = psdSetToPsfSet(PSD[:-nNaturalGS],wvl,fao.psInMas,scaleFactor=(2*np.pi*1e-9/wvl)**2)
     
        
    if doConvolve == False:
        results = psfLongExpPointingsArr
    else:
        # LOW ORDER PART
        psInMas_NGS        = fao.psInMas * (wvl_LO/wvl) #airy pattern PSF FWHM
        NGS_SR, psdArray, psfLE_NGS, NGS_FWHM_mas = psdSetToPsfSet(PSD[-nNaturalGS:],wvl_LO,psInMas_NGS,scaleFactor=(2*np.pi*1e-9/wvl_LO)**2,verbose=True)
        cartPointingCoords = np.dstack( (xxPointigs, yyPointigs) ).reshape(-1, 2)
        cartNGSCoordsList = []
        for i in range(nNaturalGS):
            cartNGSCoordsList.append(polarToCartesian(polarNGSCoords[i,:]))        
        cartNGSCoords = np.asarray(cartNGSCoordsList)
        mLO                = MavisLO(path, parametersFile, windPsdFile)
        Ctot               = mLO.computeTotalResidualMatrix(np.array(cartPointingCoords), cartNGSCoords, NGS_flux, NGS_SR, NGS_FWHM_mas)
        cov_ellipses       = mLO.ellipsesFromCovMats(Ctot)
        #print(cov_ellipses)
    
        # FINAl CONVOLUTION
        results = []
        for ellp, psfLongExp in zip(cov_ellipses, psfLongExpPointingsArr):
            results.append(convolve(psfLongExp, residualToSpectrum(ellp, wvl, N, 1/(fao.fovInPixel * fao.psInMas))))
        
    results[0].standardPlot(True)
    
    # save PSF cube in fits    
    hdul1 = fits.HDUList()
    cube =[]
    hdul1.append(fits.PrimaryHDU())
    for img in results:
        cube.append(cp.asnumpy(img.sampling))
    
    cube = np.array(cube)
    hdul1.append(fits.ImageHDU(data=cube))
    hdul1.writeto( outputDir + outputFile + '.fits',overwrite=True)
    print(cube.shape)