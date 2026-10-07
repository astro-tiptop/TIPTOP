# ----------------------------------------------------------------------------
# --- Explicit Imports to prevent Namespace Pollution ---
# ----------------------------------------------------------------------------

# Core simulation classes
from .baseSimulation import baseSimulation
from .asterismSimulation import asterismSimulation
from .asterismSimulationHo import asterismSimulationHo

# Explicitly import GPU flags from the underlying libraries
from mastsel import gpuEnabled as gpuMastsel
from p3.aoSystem import gpuEnabled as gpuP3

def gpuSelect(gpuIndex):
    """
        Select the GPU used by P3 and MASTSEL (no effect if neither runs on GPU).

        :param gpuIndex: required, GPU index; if larger than the last available index, a warning is printed and the current device is kept.
        :type gpuIndex: int
    """
    if gpuMastsel or gpuP3:
        import cupy as cp
        max_index = cp.cuda.runtime.getDeviceCount() - 1
        if gpuIndex <= max_index:
            gpu_device = cp.cuda.Device(gpuIndex)
            gpu_device.use()
        else:
            print('Trying to use GPU index ', gpuIndex, ' while max index allowed is ', max_index)
            print('Defaulting to first GPU (index 0)')


# def checkParameterFile(data2check):
#     '''
#     This function can be used to verify that the parameters in the parameter file
#     fulfill basic requirements. Currently this only support the core requirements.
#     TODO:
#         verification of the type of optionnal parameter
#         Verification of lists that need to be the same length
#         Verification of the existance of parameters whan some other parameters are set
#             For example when sensor_HO.WfsType == 'pyramid' the parameter Modulation must be set.

#     Parameters
#     ----------
#     data2check : dict
#         dictionnary containing the parameter file to be checked for requirement.

#     Returns
#     -------
#     None.

#     '''
    
#     myRequiredPar = {'telescope': {'TelescopeDiameter':8.,'Resolution': 128},
#                      'atmosphere': {'Seeing': 0.6},
#                      'sources_science': {'Wavelength': [1.6e-06], 'Zenith': [0.0],
#                                          'Azimuth': [0.0]},
#                      'sources_HO': {'Wavelength': [7.5e-07]},
#                      'sensor_science': {'PixelScale': 40, 'FieldOfView': 256},
#                      'sensor_HO': {'PixelScale': 832, 'FieldOfView': 6,
#                                    'NumberPhotons': [200.0], 'SigmaRON': 0.0},
#                      'DM': {'NumberActuators': [20], 'DmPitchs': [0.25]}}
    
#     for sec in myRequiredPar.keys():
#         if sec in data2check:
#             for opt in myRequiredPar[sec].keys():
#                 if opt in data2check[sec]:
#                     if type(data2check[sec][opt])!=type(myRequiredPar[sec][opt]):
#                         if type(myRequiredPar[sec][opt])==list:
#                             #If we want a list or an int it MUST be that
#                             raise TypeError("Parameter '{}' in section '{}' must be a type '{}'"
#                                             .format(opt,sec,type(myRequiredPar[sec][opt]).__name__))
#                         elif (type(myRequiredPar[sec][opt])==float and type(data2check[sec][opt])==list 
#                               or type(myRequiredPar[sec][opt])==int and type(data2check[sec][opt])==list):
#                             # if we want a float and there is an int instead we do not care
#                             # however it cannot be a list
#                             raise TypeError("Parameter '{}' in section '{}' must be a type '{}'"
#                                             .format(opt,sec,type(myRequiredPar[sec][opt]).__name__))
#                     if type(myRequiredPar[sec][opt])==list and not data2check[sec][opt]:
#                         raise ValueError("The list '{}' in section '{}' should not be empty"
#                                          .format(opt,sec))
#                 else:
#                     #The option is missing
#                     raise KeyError("parameter '{}' is missing from section '{}' in the parameter file"
#                                    .format(opt,sec))
#         else:
#             #the section is missing
#             raise KeyError("section '{}' is not present in the parameter file"
#                            .format(sec))

def overallSimulation(path2param, parametersFile, outputDir, outputFile, doConvolve=True,
                      doPlot=False, returnRes=False, returnMetrics=False, addSrAndFwhm=True,
                      verbose=False, getHoErrorBreakDown=False, ensquaredEnergy=False,
                      eeRadiusInMas=50, savePSDs=False, saveJson=False, gpuIndex=0,
                      exactMultiWavelengthPSD=False):
    """
        Run a full TIPTOP simulation (HO PSD, LO residuals, PSFs) from a parameter file.

        :param path2param: required, path to the folder containing the parameter file.
        :type path2param: str
        :param parametersFile: required, name of the parameter file without the extension (.ini or .yml).
        :type parametersFile: str
        :param outputDir: required, path to the folder in which to write the output.
        :type outputDir: str
        :param outputFile: required, name of the output file without the extension.
        :type outputFile: str
        :param doConvolve: optional default: True, if LO is enabled, convolve the HO PSFs with the LO residual (tip/tilt) kernel of each direction. If False, the HO PSFs are returned without this convolution.
        :type doConvolve: bool
        :param doPlot: optional default: False, display the resulting PSFs.
        :type doPlot: bool
        :param returnRes: optional default: False, return the residual errors instead of saving the results (see return).
        :type returnRes: bool
        :param returnMetrics: optional default: False, return Strehl ratio, FWHM and encircled energy within eeRadiusInMas instead of saving the results (ignored if returnRes is True).
        :type returnMetrics: bool
        :param addSrAndFwhm: optional default: True, add SR, FWHM and EE of each PSF in the header of the output fits file.
        :type addSrAndFwhm: bool
        :param verbose: optional default: False, print all messages.
        :type verbose: bool
        :param getHoErrorBreakDown: optional default: False, compute and print the HO error breakdown.
        :type getHoErrorBreakDown: bool
        :param ensquaredEnergy: optional default: False, compute the ensquared energy instead of the encircled energy.
        :type ensquaredEnergy: bool
        :param eeRadiusInMas: optional default: 50, radius used for the encircled energy (if ensquaredEnergy is True, half the side of the square).
        :type eeRadiusInMas: float
        :param savePSDs: optional default: False, also save the PSDs in the output fits file.
        :type savePSDs: bool
        :param saveJson: optional default: False, save the radial PSF profiles in ``<outputFile>1D_PSF.json``.
        :type saveJson: bool
        :param gpuIndex: optional default: 0, index of the GPU used for the simulation (if GPU is available).
        :type gpuIndex: int
        :param exactMultiWavelengthPSD: optional default: False, with more than one science wavelength, compute one exact PSD grid per wavelength, so that each PSF is identical to a single-wavelength run. If False, all wavelengths share one approximate grid. If True, memory and computation time can be significantly higher than with False: the HO PSD (including reconstructor and controller) is computed and stored once per wavelength, so the cost grows roughly with the number of wavelengths.
        :type exactMultiWavelengthPSD: bool

        :return: if returnRes, HO residual in nm RMS per science direction (and LO residual in nm RMS per direction, if LO is enabled); if returnMetrics, (sr, fwhm, ee); otherwise None, and the results are saved in ``<outputDir>/<outputFile>.fits``.
        :rtype: numpy.ndarray or tuple or None

    """

    gpuSelect(gpuIndex)

    simulation = baseSimulation(path2param, parametersFile, outputDir, outputFile, doConvolve,
                      doPlot, addSrAndFwhm, verbose, getHoErrorBreakDown, savePSDs, ensquaredEnergy,
                      eeRadiusInMas, exactMultiWavelengthPSD)
    
    simulation.doOverallSimulation()

    if saveJson:
        simulation.savePSFprofileJSON()

    if returnRes:
        if simulation.LOisOn:
            return simulation.HO_res, simulation.LO_res
        else:
            return simulation.HO_res
    elif returnMetrics:
        simulation.computeMetrics()
        return simulation.sr, simulation.fwhm, simulation.ee
    else:
        simulation.saveResults()


def asterismSelection(simulName, path2param, parametersFile, outputDir, outputFile,
                      doPlot=False, returnRes=False, returnMetrics=True, addSrAndFwhm=True,
                      verbose=False, getHoErrorBreakDown=False, ensquaredEnergy=False,
                      eeRadiusInMas=50, doConvolve=False, plotInComputeAsterisms=False,
                      progressStatus=False, gpuIndex=0):

    """
        Evaluate the asterisms defined in the ``[ASTERISM_SELECTION]`` section (requires LO).

        :param simulName: required, name of the simulation, used as prefix of the saved/reloaded .npy files.
        :type simulName: str
        :param path2param: required, path to the folder containing the parameter file.
        :type path2param: str
        :param parametersFile: required, name of the parameter file without the extension (.ini or .yml).
        :type parametersFile: str
        :param outputDir: required, path to the folder in which to write the output.
        :type outputDir: str
        :param outputFile: required, name of the output file without the extension.
        :type outputFile: str
        :param doPlot: optional default: False, display intermediate plots.
        :type doPlot: bool
        :param returnRes: optional default: False, return the HO and LO residuals (see return).
        :type returnRes: bool
        :param returnMetrics: optional default: True, return Strehl ratio, FWHM, encircled energy and covariance ellipses (ignored if returnRes is True).
        :type returnMetrics: bool
        :param addSrAndFwhm: optional default: True, add SR and FWHM in the header of the output fits file.
        :type addSrAndFwhm: bool
        :param verbose: optional default: False, print all messages.
        :type verbose: bool
        :param getHoErrorBreakDown: optional default: False, currently ignored.
        :type getHoErrorBreakDown: bool
        :param ensquaredEnergy: optional default: False, currently ignored.
        :type ensquaredEnergy: bool
        :param eeRadiusInMas: optional default: 50, radius used for the encircled energy.
        :type eeRadiusInMas: float
        :param doConvolve: optional default: False, convolve the HO PSFs with the LO residual kernel of each asterism.
        :type doConvolve: bool
        :param plotInComputeAsterisms: optional default: False, display the asterisms.
        :type plotInComputeAsterisms: bool
        :param progressStatus: optional default: False, display the progress status.
        :type progressStatus: bool
        :param gpuIndex: optional default: 0, index of the GPU used for the simulation (if GPU is available).
        :type gpuIndex: int

        :return: if returnRes, (HO_res, LO_res, simulation); if returnMetrics, (sr, fwhm, ee, cov_ellipses, simulation); otherwise simulation. None if there is no ``[ASTERISM_SELECTION]`` section or LO is not enabled.
        :rtype: tuple or asterismSimulation or None

    """

    gpuSelect(gpuIndex)

    simulation = asterismSimulation(simulName, path2param, parametersFile, outputDir, outputFile,
                      doPlot, addSrAndFwhm, verbose, progressStatus=progressStatus)


    if simulation.hasAsterismSection and simulation.LOisOn:

        simulation.computeAsterisms(eeRadiusInMas, doConvolve=doConvolve, plotGS=plotInComputeAsterisms)

        if returnRes:
            return simulation.HO_res_Asterism, simulation.LO_res_Asterism, simulation
        elif returnMetrics:
            return simulation.strehl_Asterism, simulation.fwhm_Asterism, simulation.ee_Asterism, simulation.cov_ellipses_Asterism, simulation
        else:
            return simulation
    else:
        return


def hoAsterismSelection(simulName, path2param, parametersFile, outputDir, outputFile,
                        doPlot=False, returnRes=False, returnMetrics=True, addSrAndFwhm=True,
                        verbose=False, getHoErrorBreakDown=False, ensquaredEnergy=False,
                        eeRadiusInMas=50, progressStatus=False, gpuIndex=0):
    """
        Evaluate the HO asterisms defined in the ``[HO_ASTERISM_SELECTION]`` section, using P3 only (no LO).

        :param simulName: required, name of the simulation, used as prefix of the saved/reloaded .npy files.
        :type simulName: str
        :param path2param: required, path to the folder containing the parameter file.
        :type path2param: str
        :param parametersFile: required, name of the parameter file without the extension (.ini or .yml).
        :type parametersFile: str
        :param outputDir: required, path to the folder in which to write the output.
        :type outputDir: str
        :param outputFile: required, name of the output file without the extension.
        :type outputFile: str
        :param doPlot: optional default: False, plot the results.
        :type doPlot: bool
        :param returnRes: optional default: False, return the HO residuals (see return).
        :type returnRes: bool
        :param returnMetrics: optional default: True, return Strehl ratio, FWHM and encircled energy (ignored if returnRes is True).
        :type returnMetrics: bool
        :param addSrAndFwhm: optional default: True, add SR and FWHM in the header of the output fits file.
        :type addSrAndFwhm: bool
        :param verbose: optional default: False, print all messages.
        :type verbose: bool
        :param getHoErrorBreakDown: optional default: False, compute the HO error breakdown.
        :type getHoErrorBreakDown: bool
        :param ensquaredEnergy: optional default: False, currently ignored.
        :type ensquaredEnergy: bool
        :param eeRadiusInMas: optional default: 50, radius used for the encircled energy.
        :type eeRadiusInMas: float
        :param progressStatus: optional default: False, display the progress status.
        :type progressStatus: bool
        :param gpuIndex: optional default: 0, index of the GPU used for the simulation (if GPU is available).
        :type gpuIndex: int

        :return: if returnRes, (HO_res, simulation); if returnMetrics, (sr, fwhm, ee, simulation); otherwise simulation. None if there is no ``[HO_ASTERISM_SELECTION]`` section.
        :rtype: tuple or asterismSimulationHo or None

    """

    gpuSelect(gpuIndex)
    
    simulation = asterismSimulationHo(simulName, path2param, parametersFile, outputDir, outputFile,
                                     doPlot, addSrAndFwhm, verbose, getHoErrorBreakDown, progressStatus)

    if simulation.hasHoAsterismSection:
        results = simulation.computeHoAsterisms(eeRadiusInMas)

        if doPlot:
            simulation.plotHoResults()

        if returnRes:
            return simulation.ho_res_HoAsterism, simulation
        elif returnMetrics:
            return simulation.strehl_HoAsterism, simulation.fwhm_HoAsterism, simulation.ee_HoAsterism, simulation
        else:
            return simulation
    else:
        print("No HO_ASTERISM_SELECTION section found in parameter file")
        return None


def reloadAsterismSelection(simulName, path2param, parametersFile, outputDir, outputFile,
                      doPlot=False, returnRes=False, returnMetrics=True, addSrAndFwhm=True,
                      verbose=False, getHoErrorBreakDown=False, ensquaredEnergy=False,
                      eeRadiusInMas=50, gpuIndex=0):
    """
        Reload the results of a previous asterismSelection run (the .npy files saved in outputDir with prefix simulName).

        :param simulName: required, name of the simulation, used as prefix of the saved/reloaded .npy files.
        :type simulName: str
        :param path2param: required, path to the folder containing the parameter file.
        :type path2param: str
        :param parametersFile: required, name of the parameter file without the extension (.ini or .yml).
        :type parametersFile: str
        :param outputDir: required, path to the folder in which to write the output.
        :type outputDir: str
        :param outputFile: required, name of the output file without the extension.
        :type outputFile: str
        :param doPlot: optional default: False, display plots in later processing (e.g. heuristic model fit).
        :type doPlot: bool
        :param addSrAndFwhm: optional default: True, passed to the simulation object.
        :type addSrAndFwhm: bool
        :param verbose: optional default: False, print all messages.
        :type verbose: bool
        :param getHoErrorBreakDown: optional default: False, passed to the simulation object.
        :type getHoErrorBreakDown: bool
        :param returnRes: currently ignored.
        :param returnMetrics: currently ignored.
        :param ensquaredEnergy: currently ignored.
        :param eeRadiusInMas: currently ignored.
        :param gpuIndex: optional default: 0, index of the GPU used (if GPU is available).
        :type gpuIndex: int

        :return: (sr, fwhm, ee, cov_ellipses, simulation)
        :rtype: tuple

    """

    gpuSelect(gpuIndex)

    simulation = asterismSimulation(simulName, path2param, parametersFile, outputDir, outputFile,
                                    doPlot, addSrAndFwhm, verbose, getHoErrorBreakDown)
    simulation.reloadResults()
    return simulation.strehl_Asterism, simulation.fwhm_Asterism, simulation.ee_Asterism, simulation.cov_ellipses_Asterism, simulation


def generateHeuristicModel(simulName, path2param, parametersFile, outputDir, outputFile, doPlot=False, doTest=True,
                      share = 0.9, eeRadiusInMas=50, gpuIndex=0):
    """
        Run asterismSelection, then fit a heuristic model of the asterism metrics on the first
        ``share`` fraction of the fields and, optionally, test it on the remaining ones.
        The model is saved in outputDir as ``<parametersFile>_hmodel``.

        :param simulName: required, name of the simulation, used as prefix of the saved/reloaded .npy files.
        :type simulName: str
        :param path2param: required, path to the folder containing the parameter file.
        :type path2param: str
        :param parametersFile: required, name of the parameter file without the extension (.ini or .yml).
        :type parametersFile: str
        :param outputDir: required, path to the folder in which to write the output.
        :type outputDir: str
        :param outputFile: required, name of the output file without the extension.
        :type outputFile: str
        :param doPlot: optional default: False, display the fit and test plots.
        :type doPlot: bool
        :param doTest: optional default: True, test the model on the fields not used for the fit.
        :type doTest: bool
        :param share: optional default: 0.9, fraction of the fields used for the fit.
        :type share: float
        :param eeRadiusInMas: currently ignored.
        :param gpuIndex: optional default: 0, index of the GPU used (if GPU is available).
        :type gpuIndex: int

        :return: the asterism simulation object
        :rtype: asterismSimulation

    """

    sr, fw, ee, covs, simul = asterismSelection(simulName, path2param, parametersFile, outputDir, outputFile, doPlot=False, doConvolve=False, gpuIndex=gpuIndex)

    sr, fw, ee, covs, simul = reloadAsterismSelection(simulName, path2param, parametersFile, outputDir, outputFile, doPlot=doPlot, gpuIndex=gpuIndex)

    simul.fitHeuristicModel(0, int(share*simul.nfields), parametersFile.split('.')[0]+'_hmodel')

    if doTest:
        simul.testHeuristicModel(int(share*simul.nfields), simul.nfields-1, parametersFile.split('.')[0]+'_hmodel', [])

    return simul