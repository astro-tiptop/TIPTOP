from .baseSimulation import *
from .asterismSimulation import *

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
    simulation = baseSimulation(path, parametersFile, outputDir, outputFile, doConvolve,
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
        
def asterismSelection(simulName, path, parametersFile, outputDir, outputFile,
                      doPlot=False, returnRes=False, returnMetrics=True, addSrAndFwhm=True,
                      verbose=False, getHoErrorBreakDown=False, eeRadiusInMas=50):

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
    simulation = asterismSimulation(simulName, path, parametersFile, outputDir, outputFile,
                      doPlot, addSrAndFwhm, verbose, getHoErrorBreakDown)

    
    if simulation.hasAsterismSection and simulation.LOisOn:
 
        simulation.computeAsterisms(eeRadiusInMas)
        
        if returnRes:
            return simulation.HO_res_Asterism, simulation.LO_res_Asterism, simulation
        elif returnMetrics:
            # this must be probably done inside computeAsterisms, at the end
            # simulation.computeMetrics(eeRadiusInMas)            
            return simulation.sr_Asterism, simulation.fwhm_Asterism, simulation.ee_Asterism, simulation.cov_ellipses_Asterism, simulation
        else:
            # this must be probably done inside computeAsterisms, at the end
            # simulation.saveResults()
            return simulation
    else:
        return


def reloadAsterismSelection(simulName, path, parametersFile, outputDir, outputFile,
                      doPlot=False, returnRes=False, returnMetrics=True, addSrAndFwhm=True,
                      verbose=False, getHoErrorBreakDown=False, eeRadiusInMas=50):

    simulation = asterismSimulation(simulName, path, parametersFile, outputDir, outputFile,
                                    doPlot, addSrAndFwhm, verbose, getHoErrorBreakDown)
    simulation.reloadResults()
    return simulation.sr_Asterism, simulation.fwhm_Asterism, simulation.ee_Asterism, simulation.cov_ellipses_Asterism, simulation
