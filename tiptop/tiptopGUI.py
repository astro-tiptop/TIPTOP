from mastsel.mavis import *

from ipywidgets import interact, interactive, fixed, interact_manual, Layout
import ipywidgets as widgets
from configparser import ConfigParser

sections_text = {}
sections_text['telescope'] = "Telescope"
sections_text['atmosphere'] = "Atmosphere"
sections_text['PSF_DIRECTIONS'] = "PSF Directions"
sections_text['GUIDESTARS_HO'] = "Guidestars for High Order"
sections_text['DM'] = "Deformable Mirrors"
sections_text['SENSOR_HO'] = "High Order Sensor"
sections_text['SENSOR_LO'] = "Low Order Sensor"
variables_text = {}
variables_text['TelescopeDiameter'] = "Telescope Diameter"
variables_text['zenithAngle'] = "Zenith Angle"
variables_text['obscurationRatio'] = "Obscuration Ratio"
variables_text['resolution'] = "Resolution used for pupil simulation"
variables_text['atmosphereWavelength'] = "Wavelength (atmosphere simulation)"
variables_text['seeing'] = "Seeing"
variables_text['L0'] = "L0"
variables_text['Cn2Weights'] = "Fractionl layers' weights"
variables_text['Cn2Heights'] = "Layers altitudes [m]"
variables_text['wSpeed'] = "Wind speeds [m/s]"
variables_text['wDir'] = "Wind directions [degrees]"
variables_text['nLayersReconstructed'] = "Number of reconstructed layers"
variables_text['ScienceWavelength'] = ""
variables_text['ScienceZenith'] = ""
variables_text['ScienceAzimuth'] = ""
variables_text['psInMas'] = ""
variables_text['psf_FoV'] = ""
variables_text['technical_FoV'] = ""
variables_text['GuideStarZenith_HO'] = ""
variables_text['GuideStarAzimuth_HO'] = ""
variables_text['GuideStarHeight_HO'] = ""
variables_text['DmPitchs'] = ""
variables_text['DmHeights'] = ""
variables_text['OptimizationZenith'] = ""
variables_text['OptimizationAzimuth'] = ""
variables_text['OptimizationWeight'] = ""
variables_text['OptimizationConditioning'] = ""
variables_text['nLenslet_HO'] = ""
variables_text['SensingWavelength_HO'] = ""
variables_text['loopGain_HO'] = ""
variables_text['SensorFrameRate_HO'] = ""
variables_text['loopDelaySteps_HO'] = ""
variables_text['nph_HO'] = ""
variables_text['sigmaRON_HO'] = ""
variables_text['Npix_per_subap_HO'] = ""
variables_text['pixel_scale_HO'] = ""
variables_text['spot_FWHM_HO_short'] = ""
variables_text['spot_FWHM_HO_long'] = ""
variables_text['N_sa_tot_LO'] = ""
variables_text['SensingWavelength_LO'] = ""
variables_text['SensorFrameRate_LO'] = ""
variables_text['loopDelaySteps_LO'] = ""
variables_text['pixel_scale_LO'] = ""
variables_text['Npix_per_subap_LO'] = ""
variables_text['WindowRadiusWCoG_LO'] = ""
variables_text['sigmaRON_LO'] = ""
variables_text['ExcessNoiseFactor_LO'] = ""
variables_text['Dark_LO'] = ""
variables_text['skyBackground_LO'] = ""
variables_text['ThresholdWCoG_LO'] = ""
variables_text['NewValueThrPix_LO'] = ""

variables_text['GuideStarZenith_LO'] = ""
variables_text['GuideStarAzimuth_LO'] = ""

variables_text['nph_LO'] = ""

variables_text['nDms'] = ""
variables_text['nOpt'] = ""
variables_text['Cn2RefHeight'] = ""

variables_text['path_pupil'] = ""
variables_text['path_static'] = ""
variables_text['oneWindSpeed'] = ""
variables_text['r0_Value'] = ""


variables_min_max_type = {}
variables_min_max_type['TelescopeDiameter'] = (0, 50, 1, float)
variables_min_max_type['zenithAngle'] = (0, 360, 10, float)
variables_min_max_type['obscurationRatio'] = (0, 1, 0.1, float)
variables_min_max_type['resolution'] = (32, 2048, 100, int)
variables_min_max_type['atmosphereWavelength'] = (1e-9, 1e-5,  100e-9, float)
variables_min_max_type['seeing'] = (0,1, 0.1, float)
variables_min_max_type['L0'] = (0,1000, 1, float)
variables_min_max_type['Cn2Weights'] = (0,1, 0.1, float)
variables_min_max_type['Cn2Heights'] = (0,100000, 1000, float)
variables_min_max_type['Cn2RefHeight'] = (0,100000, 1000, float)
variables_min_max_type['wSpeed'] = (0, 1000, 1, float)
variables_min_max_type['wDir'] = (0, 360, 10, float)
variables_min_max_type['nLayersReconstructed'] = (1, 100, 1, int)
variables_min_max_type['ScienceWavelength'] = (1e-9, 1e-5, 100e-9, float)
variables_min_max_type['ScienceZenith'] = (0, 360, 10, float)
variables_min_max_type['ScienceAzimuth'] = (0, 360, 10,  float)
variables_min_max_type['psInMas'] = (0,1000, 1, float)
variables_min_max_type['psf_FoV'] = (0, 1000, 1, float)
variables_min_max_type['technical_FoV'] = (1, 1000, 10, float)
variables_min_max_type['GuideStarZenith_HO'] = (0, 360, 10, float)
variables_min_max_type['GuideStarAzimuth_HO'] = (0, 360, 10, float)
variables_min_max_type['GuideStarHeight_HO'] = (0, 500000, 1000, float)
variables_min_max_type['DmPitchs'] = (0, 1, 0.1, float)
variables_min_max_type['DmHeights'] = (0, 100000, 1000, float)
variables_min_max_type['OptimizationZenith'] = (0, 360, 10, float)
variables_min_max_type['OptimizationAzimuth'] = (0, 360, 10, float)
variables_min_max_type['OptimizationWeight'] = (0,100, 1, float)
variables_min_max_type['OptimizationConditioning'] = (0, 1000, 1, float)
variables_min_max_type['nLenslet_HO'] = (1, 1000, 1, int)
variables_min_max_type['SensingWavelength_HO'] = (1e-9, 1e-5, 100e-9, float)
variables_min_max_type['loopGain_HO'] = (0, 1, 0.1, float)
variables_min_max_type['SensorFrameRate_HO'] = (1, 10000, 100, float)
variables_min_max_type['loopDelaySteps_HO'] = (0, 100, 1, int)
variables_min_max_type['nph_HO'] = (0, 10000, 1, float)
variables_min_max_type['sigmaRON_HO'] = (0, 10000, 0.1, float)
variables_min_max_type['Npix_per_subap_HO'] = (1, 10000, 1, int)
variables_min_max_type['pixel_scale_HO'] = (0.1, 10000, 1, float)
variables_min_max_type['spot_FWHM_HO_short'] = (1, 10000, 100, int)
variables_min_max_type['spot_FWHM_HO_long'] = (1, 10000, 100, int)
variables_min_max_type['N_sa_tot_LO'] = (1, 100, 1, int)
variables_min_max_type['SensingWavelength_LO'] = (1e-9, 1e-5, 100e-9, float)
variables_min_max_type['SensorFrameRate_LO'] = (1, 10000, 100, float)
variables_min_max_type['loopDelaySteps_LO'] = (0, 100, 1, int)
variables_min_max_type['pixel_scale_LO'] = (1, 1000, 1, float)
variables_min_max_type['Npix_per_subap_LO'] = (1, 1000, 1, int)
variables_min_max_type['WindowRadiusWCoG_LO'] = (1, 10, 1, int)
variables_min_max_type['sigmaRON_LO'] = (0, 1000,  0.1, float)
variables_min_max_type['ExcessNoiseFactor_LO'] = (0, 1000, 0.1, float)
variables_min_max_type['Dark_LO'] = (0, 1000, 0.1, float)
variables_min_max_type['skyBackground_LO'] = (0, 1000, 0.1, float)
variables_min_max_type['ThresholdWCoG_LO'] = (0, 1000, 0.1, float)
variables_min_max_type['NewValueThrPix_LO'] = (0, 1000, 0.1, float)

variables_min_max_type['GuideStarZenith_LO'] = (0, 360, 1.0, float)
variables_min_max_type['GuideStarAzimuth_LO'] = (0, 360, 1.0, float)
variables_min_max_type['nph_LO'] = (0, 200, 1, int)

variables_min_max_type['nDms'] = (1, 30, 1, int)
variables_min_max_type['nOpt'] = (1, 30, 1, int)

variables_min_max_type['path_pupil'] = ("", "", "", str)
variables_min_max_type['path_static'] = ("", "", "", str)
variables_min_max_type['oneWindSpeed'] = (0, 100, 1.0, float)
variables_min_max_type['r0_Value'] = (0, 1000, 1.0, float)
