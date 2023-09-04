import os
import numpy as np
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

def psdSetToPsfSet(N, freq_range, dk, mask, inputPSDs,wavelength,pixelscale,npixel,scaleFactor=1,verbose=False):
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
        if verbose:
            print('SR(@',int(wavelength*1e9),'nm)  :', SR)
            print('FWHM(@',int(wavelength*1e9),'nm):', FWHM)

    return NGS_SR, psdArray, psfLongExpArr, NGS_FWHM_mas

def overallSimulation(path, parametersFile, outputDir, outputFile, doConvolve=False,
                      doPlot=False, returnRes=False, addSrAndFwhm=False,
                      verbose=False, getHoErrorBreakDown=False, savePSDs=False):
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
    
    if returnRes and doPlot:
        print('WARNING: returnRes and doPlot cannot both be True, setting doPlot to False.')
        doPlot = False
        
    # initiate the parser
    fullPathFilename_ini = os.path.join(path, parametersFile + '.ini')
    fullPathFilename_yml = os.path.join(path, parametersFile + '.yml')
    
    # initialize jitter_FWHM variable with a default value
    jitter_FWHM = None

    if os.path.exists(fullPathFilename_yml):
        fullPathFilename = fullPathFilename_yml
        with open(fullPathFilename_yml) as f:
            my_yaml_dict = yaml.safe_load(f)        
        my_data_map = my_yaml_dict
    elif os.path.exists(fullPathFilename_ini):
        fullPathFilename = fullPathFilename_ini
        config = ConfigParser()
        config.optionxform = str
        config.read(fullPathFilename_ini)
        my_data_map = {} 
        for section in config.sections():
            my_data_map[section] = {}
            for name,value in config.items(section):
                my_data_map[section].update({name:eval(value)})
    else:
        raise FileNotFoundError('No .yml or .ini can be found in '+ path)
                  
    # read main parameters
    tel_radius = my_data_map['telescope']['TelescopeDiameter']/2  # mas
    wvl_temp = my_data_map['sources_science']['Wavelength']
    if isinstance(wvl_temp, list):
        wvl = wvl_temp[0]  # lambda
    else:
        wvl = wvl_temp     # lambda
    zenithSrc  = my_data_map['sources_science']['Zenith']
    azimuthSrc = my_data_map['sources_science']['Azimuth']

    # it checks if LO parameters are set and then it acts accordingly
    if 'sensor_LO' in my_data_map.keys():
        LOisOn = True
        if verbose: print('LO part is present')
    else:
        LOisOn = False
        nNaturalGS = 0
        if verbose: print('LO part is not present')

    if LOisOn:
        LO_wvl_temp = my_data_map['sources_LO']['Wavelength']
        if isinstance(LO_wvl_temp, list):
            LO_wvl = LO_wvl_temp[0]  # lambda
        else:
            LO_wvl = LO_wvl_temp     # lambda
        LO_zen     = my_data_map['sources_LO']['Zenith']
        LO_az      = my_data_map['sources_LO']['Azimuth']
        LO_fluxes  = my_data_map['sensor_LO']['NumberPhotons']
        fr         = my_data_map['RTC']['SensorFrameRate_LO']

    if 'jitter_FWHM' in my_data_map['telescope'].keys():
        jitter_FWHM = my_data_map['telescope']['jitter_FWHM']

    fao = fourierModel( fullPathFilename, calcPSF=False, verbose=verbose
                       , display=False, getPSDatNGSpositions=True
                       , computeFocalAnisoCov=False, TiltFilter=LOisOn
                       , getErrorBreakDown=getHoErrorBreakDown)

    if LOisOn:
        # NGSs positions
        NGS_flux = []
        polarNGSCoordsList = []
        for aFlux, aZen, aAz in zip(LO_fluxes, LO_zen, LO_az):
            polarNGSCoordsList.append([aZen, aAz])
            NGS_flux.append(aFlux*fr)
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
    nPixPSF       = int(fao.freq.nOtf /fao.freq.kRef_)
    freq_range    = fao.ao.cam.fovInPix*fao.freq.PSDstep
    pitch         = 1/freq_range
    grid_diameter = pitch*N
    sx            = int(2*np.round(tel_radius/pitch))
    dk            = 1e9*fao.freq.kcMax_/fao.freq.resAO

    # Define the pupil shape
    mask = Field(wvl, N, grid_diameter)

    psInMas = cpuArray(fao.freq.psInMas)

    mask.sampling = congrid(arrayP3toMastsel(fao.ao.tel.pupil), [sx, sx])
        
    mask.sampling = zeroPad(mask.sampling, (N-sx)//2)
    
    if verbose:
        print('fao.samp:', fao.freq.samp)
        print('fao.PSD.shape:', fao.PSD.shape)
        print('fao.freq.psInMas:', psInMas)

    # HO PSF
    if verbose:
        print('******** HO PSF')
    
    pointings_SR, psdPointingsArray, psfLongExpPointingsArr, pointings_FWHM_mas = psdSetToPsfSet(N, 
                                                                                                 freq_range, 
                                                                                                 dk,
                                                                                                 mask, 
                                                                                                 arrayP3toMastsel(PSD[0:nPointings]),
                                                                                                 wvl,
                                                                                                 psInMas[0],
                                                                                                 nPixPSF,
                                                                                                 scaleFactor=(2*np.pi*1e-9/wvl)**2,
                                                                                                 verbose=verbose)

    if LOisOn == False or (doConvolve == False and returnRes == False):
        results = []
        for psfLongExp in psfLongExpPointingsArr:
            if jitter_FWHM is not None:
                ellp = [0, sigma_from_FWHM(jitter_FWHM), sigma_from_FWHM(jitter_FWHM)]
                results.append(convolve(psfLongExp,
                               residualToSpectrum(ellp, wvl, N, 1/(fao.ao.cam.fovInPix * psInMas[0]))))
            else:
                results.append(psfLongExp)

    else:
        # LOW ORDER PART
        if verbose:
            print('******** LO PART')
        psInMas_NGS        = psInMas[0] * (LO_wvl/wvl) #airy pattern PSF FWHM

        if verbose:
            print('******** HO PSF - NGS directions')
        NGS_SR, psdArray, psfLE_NGS, NGS_FWHM_mas = psdSetToPsfSet(N, 
                                                                   freq_range, 
                                                                   dk, 
                                                                   mask, 
                                                                   arrayP3toMastsel(PSD[-nNaturalGS:]), 
                                                                   LO_wvl, 
                                                                   psInMas_NGS, 
                                                                   nPixPSF,
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
                resSpec = residualToSpectrum(ellp, wvl, N, 1/(fao.ao.cam.fovInPix * psInMas[0]))
                results.append(convolve(psfLongExp, resSpec))

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
        psdOL.sampling = arrayP3toMastsel(fao.ao.atm.spectrum(k) * pf * (fao.freq.wvlRef/np.pi)**2) # the PSD must be provided in m^2.m^2
        # Get the OPEN-LOOP PSF
        psfOL = longExposurePsf(mask, psdOL)
        # It cuts the PSF if the PSF is larger than the requested dimension (N>nPixPSF)
#        if psfOL.sampling.shape[0] > nPixPSF:
#            psfOL.sampling = psfOL.sampling[int(psfOL.sampling.shape[0]/2-nPixPSF/2):int(psfOL.sampling.shape[0]/2+nPixPSF/2),
#                                            int(psfOL.sampling.shape[1]/2-nPixPSF/2):int(psfOL.sampling.shape[1]/2+nPixPSF/2)]
#        psfOL.sampling = cpuArray(psdOL.sampling)

        if doPlot:
            fig, ax1 = plt.subplots(1,1)
            im = ax1.imshow(np.log(np.abs(psfOL.sampling) + 1e-20), cmap='hot')
            ax1.set_title('open loop PSF', color='black')
            
        # DIFFRACTION LIMITED PSD
        psdDL = Field(wvl, N, freq_range, 'rad')
        psfDL = longExposurePsf(mask, psdDL)
        # It cuts the PSF if the PSF is larger than the requested dimension (N>nPixPSF)
        if psfDL.sampling.shape[0] > nPixPSF:
            psfDL.sampling = psfDL.sampling[int(psfDL.sampling.shape[0]/2-nPixPSF/2):int(psfDL.sampling.shape[0]/2+nPixPSF/2),
                                            int(psfDL.sampling.shape[1]/2-nPixPSF/2):int(psfDL.sampling.shape[1]/2+nPixPSF/2)]
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
        if savePSDs:
            hdul1.append(fits.ImageHDU(data=cpuArray(PSD))) # append high order PSD


        #############################
        # header
        hdr0 = hdul1[0].header
        now = datetime.now()
        hdr0['TIME'] = now.strftime("%Y%m%d_%H%M%S")
        # parameters in the header
        for key_primary in my_data_map:
            for key_secondary in my_data_map[key_primary]:
                temp = my_data_map[key_primary][key_secondary]
                if isinstance(temp, list):
                    iii = 0
                    for elem in temp:
                        if isinstance(elem, list):
                            jjj = 0
                            for elem2 in elem:
                                hdr0['HIERARCH '+key_primary+' '+key_secondary +' '+str(iii)+' '+str(jjj)] = elem2
                                jjj += 1
                        else:                        
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
        hdr1['PIX_MAS'] = str(psInMas[0])
        hdr1['CC'] = "CARTESIAN COORD. IN ASEC OF THE "+str(pp.shape[1])+" SOURCES"
        for i in range(pp.shape[1]):
            hdr1['CCX'+str(i).zfill(4)] = pp[0,i]
            hdr1['CCY'+str(i).zfill(4)] = pp[1,i]
        if addSrAndFwhm:
            for i in range(cube.shape[0]):
                hdr1['SR'+str(i).zfill(4)]   = getStrehl(cube[i,:,:], fao.ao.tel.pupil, fao.freq.sampRef)
            for i in range(cube.shape[0]):
                hdr1['FWHM'+str(i).zfill(4)] = getFWHM(cube[i,:,:], psInMas[0], method='contour', nargout=1)

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

        if savePSDs:
            # header of the PSD
            hdr4 = hdul1[4].header
            hdr4['TIME'] = now.strftime("%Y%m%d_%H%M%S")
            hdr4['CONTENT'] = "High Order PSD"
            hdr4['SIZE'] = str(PSD.shape)

        #############################

        hdul1.writeto( os.path.join(outputDir, outputFile + '.fits'), overwrite=True)
        if verbose:
            print("Output cube shape:", cube.shape)
