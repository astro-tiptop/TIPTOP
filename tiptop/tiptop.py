from .baseSimulation import *
from .asterismSimulation import *

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
                      doPlot=False, returnRes=False, returnMetrics=False, addSrAndFwhm=False,
                      verbose=False, getHoErrorBreakDown=False, ensquaredEnergy=False,
                      eeRadiusInMas=50, savePSDs=False):
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
        :param ensquaredEnergy: optional default: False, If you want ensquared energy instead of encircled energy set this to True.
        :type ensquaredEnergy: bool
        :param eeRadiusInMas: optional default: 50, used together with returnMetrics, radius used for the computation of the encirlced energy (if ensquaredEnergy is selected, this is half the side of the square)
        :type eeRadiusInMas: float
        :param savePSDs: optional default: False, If you want to save PSD in the output fits file set this to True.
        :type savePSDs: bool

        :return: TBD
        :rtype: TBD

    """
    simulation = baseSimulation(path2param, parametersFile, outputDir, outputFile, doConvolve,
                      doPlot, addSrAndFwhm, verbose, getHoErrorBreakDown, savePSDs, ensquaredEnergy,
                      eeRadiusInMas)
    
    simulation.doOverallSimulation()

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
                      progressStatus=False):

    """
        function to run the entire tiptop asterism evaluation on the input file

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
        :param ensquaredEnergy: optional default: False, If you want ensquared energy instead of encircled energy set this to True.
        :type ensquaredEnergy: bool
        :param eeRadiusInMas: optional default: 50, used together with returnMetrics, radius used for the computation of the encirlced energy
        :type eeRadiusInMas: float
        :param plotInComputeAsterisms: optional default: False, If you want to display asterisms.
        :type plotInComputeAsterisms: bool
        :param progressStatus: optional default: False, If you want to display progress status.
        :type progressStatus: bool

        :return: TBD
        :rtype: TBD

    """
    simulation = asterismSimulation(simulName, path2param, parametersFile, outputDir, outputFile,
                      doPlot, addSrAndFwhm, verbose, progressStatus=progressStatus)

    
    if simulation.hasAsterismSection and simulation.LOisOn:
 
        simulation.computeAsterisms(eeRadiusInMas, doConvolve=doConvolve, plotGS=plotInComputeAsterisms)
        
        if returnRes:
            return simulation.HO_res_Asterism, simulation.LO_res_Asterism, simulation
        elif returnMetrics:
            # this must be probably done inside computeAsterisms, at the end
            # simulation.computeMetrics(eeRadiusInMas)            
            return simulation.strehl_Asterism, simulation.fwhm_Asterism, simulation.ee_Asterism, simulation.cov_ellipses_Asterism, simulation
        else:
            # this must be probably done inside computeAsterisms, at the end
            # simulation.saveResults()
            return simulation
    else:
        return


def reloadAsterismSelection(simulName, path2param, parametersFile, outputDir, outputFile,
                      doPlot=False, returnRes=False, returnMetrics=True, addSrAndFwhm=True,
                      verbose=False, getHoErrorBreakDown=False, ensquaredEnergy=False,
                      eeRadiusInMas=50):
    simulation = asterismSimulation(simulName, path2param, parametersFile, outputDir, outputFile,
                                    doPlot, addSrAndFwhm, verbose, getHoErrorBreakDown)
    simulation.reloadResults()
    return simulation.strehl_Asterism, simulation.fwhm_Asterism, simulation.ee_Asterism, simulation.cov_ellipses_Asterism, simulation
