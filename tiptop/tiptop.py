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
