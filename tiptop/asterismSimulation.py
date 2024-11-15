import matplotlib
from string import ascii_uppercase
from .baseSimulation import *
import math
from scipy.optimize import curve_fit, minimize
from scipy import interpolate

from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import LogStretch

from tiptop.nnModel import *
import mpl_scatter_density

import pickle
import os.path
import sys
import time

from dataclasses import dataclass
from typing import List

@dataclass
class Star:
    zenith: float
    azimuth: float
    photons: float
    freq: float

@dataclass
class AsterismProperties:
    index: int
    asterism: List[Star]
    jitter: float
    strehl: float
    fwhm: float
    encircled_energy: float

def funcPolar(X, A0, A, B, C, D, E0, E, F, G, H, I, J, J0, J1):
    r, af = X
    f = af/100
    return (A0 * (r/100) ** 3) + (A * (r/100) ** 2) + (B * r/100 ) + (C * r ** 0.5) + (D * r ** 0.25) + (E0 * (0.1/f) ** 3) + (E * (0.1/f) ** 2) + (F * 1/f ) + (G * (100/f) ** 0.5) + (H * (100/f) ** 0.25) + I + (J * (r/f)**2) + (J0 * (r/f)) + (J1 * (r/f)**0.5) 

def costFunction(params, x, y):
    A0, A, B, C, D, E, F, G, H, I, J = params
    return np.sum( ( (y - funcPolar(x, A0, A, B, C, D, E, F, G, H, I, J)) / y ) **2 )


def funcCartesian(X, A, B, C, D, E, F, G, H, I, J, K, L):
    r, x, y, f = X
    #f = np.log(1/af)
    return (A * x ** 3) + (B * y ** 3) + (C * f ** 3) + \
           (D * x ** 2) + (E * y ** 2) + (F * f ** 2) + \
           (G * x ) + (H * y ) + (I * f ) + \
           (J * x * y) + K + L*(x * y * f)

def funcMix(X, A, B, C, D, E, F, G, H, I, J):
    r, x, y, f = X
    #f = np.log(1/af)
    return (A * r ** 3) + (B * f ** 3) + \
           (C * r ** 2) + (D * f ** 2) + \
           (E * r) + (F *f) + \
           (G * y/x) + (H * (x/y)) + (I * r * f) + J

def unrollAsterismData(all_combos, c1, c2, flux, freq):
    asterism = np.array( [np.take(c1, all_combos),
                        np.take(c2, all_combos),
                        np.take(flux, all_combos),
                        np.take(freq, all_combos)] )
    return np.swapaxes(asterism, 0,1)

class asterismSimulation(baseSimulation):


    def __init__(self, simulName, path, parametersFile, outputDir,
                 outputFile, doPlot=False, addSrAndFwhm=False, verbose=False,
                 getHoErrorBreakDown=False, progressStatus=False):
        super().__init__(path, parametersFile, outputDir, outputFile, doConvolve=True,
                          doPlot=False, addSrAndFwhm=addSrAndFwhm,
                          verbose=verbose, getHoErrorBreakDown=getHoErrorBreakDown,
                          savePSDs=False)
        self.nNGS = 0
        self.firstConfigCall = True
        self.simulName = simulName
        self.doPlotAst = doPlot
        # store the parameters in data members
        self.asterismsInputDataCartesian = None
        self.asterismsInputDataPolar = None
        self.hasAsterismSection = False
        self.progressStatus = progressStatus
        if 'ASTERISM_SELECTION' in self.my_data_map.keys():
            # some global settings which are used when generating random data
            self.hasAsterismSection = True
            self.globalOffset = 0
            self.asterismMode = self.my_data_map['ASTERISM_SELECTION']['mode']
            listZ = self.my_data_map['ASTERISM_SELECTION']['Zenith']
            listA = self.my_data_map['ASTERISM_SELECTION']['Azimuth']
            listP = self.my_data_map['ASTERISM_SELECTION']['NumberPhotons']
            listF = self.my_data_map['ASTERISM_SELECTION']['Frequencies']
            self.transmissionFactor = self.my_data_map['ASTERISM_SELECTION']['transmissionFactor']
            self.bands              = self.my_data_map['ASTERISM_SELECTION']['bands']
            self.ObscurationRatio   = self.my_data_map['telescope']['ObscurationRatio']
            self.TelescopeDiameter  = self.my_data_map['telescope']['TelescopeDiameter']
            self.NumberLenslets     = self.my_data_map['sensor_LO']['NumberLenslets']
            self.N_sa_tot_LO        = self.NumberLenslets[0]**2
            if self.NumberLenslets[0] > 2:
                self.N_sa_tot_LO   = int ( np.floor( self.N_sa_tot_LO * np.pi/4.0 * (1.0 - self.ObscurationRatio**2) ) )
            self.fluxScaling = ((self.TelescopeDiameter/2.0)**2 * np.pi * (1-self.ObscurationRatio**2) * self.transmissionFactor / self.N_sa_tot_LO)
            if not isinstance(listF, list):
                listF  = [listF] * len(listZ)
            self.cumAstSizes = [0]
            self.cumStarSizes = [0]
            self.nfields = 1
            if self.asterismMode=='Sets':
                self.generateFromList(listZ, listA, listP, listF)
            elif self.asterismMode[:7]=='Singles':
                listF = [250] * 13
                if self.asterismMode[7]=='3' or self.asterismMode[7]=='1':
                    nStars = len(self.my_data_map['ASTERISM_SELECTION']['Zenith'])
                    self.nNGS = int(self.asterismMode[7])                    
                    pointings = polarToCartesian(np.array( [listZ, listA]))
                    xxPointigs  = pointings[0,:]
                    yyPointigs  = pointings[1,:]
                    self.isMono = self.nNGS==1
                else:
                    self.asterismMode = 'INVALID'
                    self.hasAsterismSection = False
                    print('ERROR: Only Singles1 (One Start Asterisms) and Singles3 (3 Stars Asterisms) are implemented.')
                    return
                all_combos = list(itertools.combinations(list(range(nStars)), self.nNGS))
                self.nfieldsSizes = [len(all_combos)]
                self.cumStarSizes.append(nStars)
                self.cumAstSizes.append(self.nfieldsSizes[0])
                zenith = np.array(listZ, dtype=np.float64)
                azimuth = np.array(listA, dtype=np.float64)
                flux = np.array(listP, dtype=np.float64)
                freq = np.array(listF, dtype=np.float64)
                self.asterismsInputDataCartesian = unrollAsterismData(all_combos, xxPointigs, yyPointigs, flux, freq)
                self.asterismsInputDataPolar = unrollAsterismData(all_combos, zenith, azimuth, flux, freq)
                self.allAsterismsIndices = self.currentFieldAsterismsIndices = np.asarray(all_combos)
            elif self.asterismMode=='Generate':
                self.asterismsInputDataCartesian, self.asterismsInputDataPolar = self.generateTriangles(listZ[0], listZ[0], 10)
            elif self.asterismMode[:4]=='File':
                self.nfields = self.my_data_map['ASTERISM_SELECTION']['fieldsNumber']
                # number of asterisms for each field
                self.nfieldsSizes = []
                self.file_field_simul = self.my_data_map['ASTERISM_SELECTION']['filename']
                self.globalOffset = self.my_data_map['ASTERISM_SELECTION']['offset']
                self.isMono = self.asterismMode[-4:]=='Mono' or self.nNGS==1
                if self.isMono:
                    self.magnitudesRange = [11,20]
                    self.fovRange = [30,30]
                    self.minStars = 1
                    self.maxStars = 20
                else:
                    self.magnitudesRange = [11,20]
                    self.fovRange = [60,60]
                    self.minStars = 3
                    self.maxStars = 12
                if self.asterismMode[4:10]=='Random':
                    self.generate_data = False
                    datafiles = []
                    datafiles.append(os.path.join(self.outputDir, self.file_field_simul +'C.npy'))
                    datafiles.append(os.path.join(self.outputDir, self.file_field_simul +'P.npy'))
                    datafiles.append(os.path.join(self.outputDir, self.file_field_simul +'F.npy'))
                    datafiles.append(os.path.join(self.outputDir, self.file_field_simul +'S.npy'))
                    datafiles.append(os.path.join(self.outputDir, self.file_field_simul +'ST.npy'))
                    datafiles.append(os.path.join(self.outputDir, self.file_field_simul +'IDX.npy'))
                    for datafile in datafiles:
                        if not os.path.exists(datafile):
                            self.generate_data = True
                    if self.generate_data:
                        self.generateRandom(self.nfields)
                    else:
                        self.asterismsInputDataCartesian = np.load(datafiles[0])
                        self.asterismsInputDataPolar = np.load(datafiles[1])
                        self.nfieldsSizes = np.load(datafiles[2]).tolist()
                        self.cumAstSizes = np.load(datafiles[3]).tolist()
                        self.cumStarSizes = np.load(datafiles[4]).tolist()
                        self.allAsterismsIndices = np.load(datafiles[5])                        
                else:
                    field_simul_data = np.load(self.file_field_simul, allow_pickle=True)
                    self.asterismsRecArray = field_simul_data
                    print('Number of Fields:', len(self.asterismsRecArray))
                    if self.isMono:
                        self.generateFromRecArray(self.nfields)
                    else:
                        self.generateFromRecArrayMulti(self.nfields)
            else:
                self.asterismMode = 'INVALID'                
        else:
            self.hasAsterismSection = False

        if 'heuristicModel' in self.my_data_map['ASTERISM_SELECTION'].keys():
            self.heuristicModel = self.my_data_map['ASTERISM_SELECTION']['heuristicModel']
            if self.isMono:
                print(os.path.join(self.outputDir, self.heuristicModel +'.npy'))
                self.monoModel = np.load(os.path.join(self.outputDir, self.heuristicModel +'.npy'), allow_pickle=True)
            # else: 
            #     load the NN here
        else:
            self.heuristicModel = None

        if self.verbose:
            print('self.cumAstSizes', self.cumAstSizes)

    def resetFieldsData(self):
        self.setsList = []
        self.polarSetsList = []
        self.nfieldsSizes = []
        self.cumAstSizes = [0]
        self.cumStarSizes = [0]
        self.skippedFieldIndexes = []
        self.allAsterismsIndices = []

    def reset_currentFieldsSourcesData(self):
        self.currentFieldsSourcesData = {}
        self.currentFieldsSourcesData['Zenith'] = []
        self.currentFieldsSourcesData['Azimuth'] = []
        self.currentFieldsSourcesData['NumberPhotons'] = []
        self.currentFieldsSourcesData['Frequencies'] = []
        self.currentFieldAsterismsIndices = []

    def appendSource(self, source):
        self.currentFieldsSourcesData['Zenith'].append(source[0])
        self.currentFieldsSourcesData['Azimuth'].append(source[1])
        self.currentFieldsSourcesData['NumberPhotons'].append(source[2])
        self.currentFieldsSourcesData['Frequencies'].append(source[3])

    def updateAsterismIndices(self, all_combos, number_of_asterisms, number_of_stars):
        self.allAsterismsIndices.extend(all_combos)
        self.nfieldsSizes.append(number_of_asterisms)
        self.cumAstSizes.append(self.cumAstSizes[-1]+number_of_asterisms)
        self.cumStarSizes.append(self.cumStarSizes[-1]+number_of_stars)

    def addFieldDataCombos(self, all_combos, number_of_asterisms, number_of_stars):
       
        self.updateAsterismIndices([*all_combos], number_of_asterisms, number_of_stars)
        
        pcoord_r = self.currentFieldsSourcesData['Zenith']
        pcoord_a = self.currentFieldsSourcesData['Azimuth']
        pcoords = np.array([pcoord_r, pcoord_a])
        rcoords = polarToCartesian(pcoords)
        xxPointigs = rcoords[0, :]
        yyPointigs = rcoords[1, :]

        flux = self.currentFieldsSourcesData['NumberPhotons']
        freq = self.currentFieldsSourcesData['Frequencies']
        
        asterism = unrollAsterismData(all_combos, xxPointigs, yyPointigs, flux, freq)
        pasterism = unrollAsterismData(all_combos, pcoords[0,:], pcoords[1,:], flux, freq)

        self.setsList.extend([*asterism])
        self.polarSetsList.extend([*pasterism])

    def generateFromList(self, listZ, listA, listP, listF):
        self.resetFieldsData()
        self.reset_currentFieldsSourcesData()
        coords = polarToCartesian(np.array( [listZ, listA]))
        xcoords  = coords[0,:]
        ycoords  = coords[1,:]
        fluxes = np.array( [listP])[0]
        freqs = np.array( [listF])[0]
        self.nNGS = coords.shape[2]
        self.currentFieldsize = coords.shape[1]
        self.nfieldsSizes = [coords.shape[1]]
        number_of_asterisms = self.currentFieldsize
        number_of_stars = 0
        all_combos = []
        setsList = []
        polarSetsList = []
        for j in range(number_of_asterisms):
            s_index = []
            for si in range(3):
                pcoords = cartesianToPolar( np.asarray([xcoords[j][si], ycoords[j][si]]))
                source = np.array([pcoords[0], pcoords[1], fluxes[j][si], freqs[j][si]])
                ss = self.sourceIsPresent(source)
                if ss==-1:
                    self.appendSource(source)
                    ss = number_of_stars
                    number_of_stars += 1
                s_index.append(ss)
            all_combos.extend([s_index])
            pcoords = cartesianToPolar( np.asarray([xcoords[j], ycoords[j]]) )
            asterism = np.vstack( [xcoords[j], ycoords[j], fluxes[j], freqs[j]] )
            pasterism = np.vstack( [pcoords[0,:], pcoords[1,:], fluxes[j], freqs[j]] )
            setsList.append(asterism)
            polarSetsList.append(pasterism)
        self.addFieldDataCombos(all_combos, number_of_asterisms, number_of_stars)
        self.asterismsInputDataCartesian = np.array(self.setsList)
        self.asterismsInputDataPolar = np.array(self.polarSetsList)
        self.allAsterismsIndices = np.array(self.allAsterismsIndices)

    def configLO(self, astIndex=None):
        if self.firstConfigCall:
            self.cartSciencePointingCoords = np.dstack( (self.xxSciencePointigs, self.yySciencePointigs) ).reshape(-1, 2)
            # Here we assume the same wavelenght for all the phon counts of the stars in the asterism
            LO_wvl_temp = self.my_data_map['sources_LO']['Wavelength']
            if isinstance(LO_wvl_temp, list):
                self.LO_wvl = LO_wvl_temp[0]  # lambda
            else:
                self.LO_wvl = LO_wvl_temp     # lambda
                
            self.my_data_map['sources_LO']['Zenith'] = self.currentFieldsSourcesData['Zenith']
            self.my_data_map['sources_LO']['Azimuth'] = self.currentFieldsSourcesData['Azimuth']
            self.my_data_map['sensor_LO']['NumberPhotons'] = self.currentFieldsSourcesData['NumberPhotons']
            self.my_data_map['RTC']['SensorFrameRate_LO'] = self.currentFieldsSourcesData['Frequencies']
            
            self.LO_psInMas      = self.my_data_map['sensor_LO']['PixelScale']
            self.LO_zen_field    = self.my_data_map['sources_LO']['Zenith']
            self.LO_az_field     = self.my_data_map['sources_LO']['Azimuth']
            self.LO_fluxes_field = self.my_data_map['sensor_LO']['NumberPhotons']
            self.LO_freqs_field  = self.my_data_map['RTC']['SensorFrameRate_LO']

            if self.check_section_key('sensor_Focus'):
                self.Focus_fluxes4s_field = self.my_data_map['sensor_Focus']['NumberPhotons']
                self.Focus_psInMas         = self.my_data_map['sensor_Focus']['PixelScale']
                Focus_wvl_temp = self.my_data_map['sources_LO']['Wavelength']
                if isinstance(Focus_wvl_temp, list):
                    self.Focus_wvl = Focus_wvl_temp[0]  # lambda
                else:
                    self.Focus_wvl = Focus_wvl_temp     # lambda
            else:
                self.Focus_fluxes4s_field = self.LO_fluxes_field
                self.Focus_psInMas         = self.LO_psInMas
                self.Focus_wvl             = self.LO_wvl
            if self.check_config_key('RTC','SensorFrameRate_Focus'):
                self.Focus_freqs_field  = self.my_data_map['RTC']['SensorFrameRate_Focus']
                if not isinstance(self.Focus_freqs_field, list):
                    self.Focus_freqs_field  = [self.Focus_freqs_field] * len(self.LO_zen_field)
            else:
                self.Focus_freqs_field = self.LO_freqs_field
        
            self.NGS_fluxes_field = []
            polarNGSCoordsList = []
            for aFr, aFlux, aZen, aAz in zip(self.LO_freqs_field, self.LO_fluxes_field, self.LO_zen_field, self.LO_az_field):
                polarNGSCoordsList.append([aZen, aAz])
                self.NGS_fluxes_field.append(aFlux*aFr)
            self.Focus_fluxes_field = []
            for aFrF, aFluxF in zip(self.Focus_freqs_field, self.Focus_fluxes4s_field):
                self.Focus_fluxes_field.append(aFluxF*aFrF)
            polarNGSCoords     = np.asarray(polarNGSCoordsList)
            self.nNaturalGS_field  = len(self.LO_zen_field)
            cartNGSCoordsList = []
            for i in range(self.nNaturalGS_field):
                cartNGSCoordsList.append(polarToCartesian(polarNGSCoords[i,:]))
            self.cartNGSCoords_field = np.asarray(cartNGSCoordsList)
        self.currentAsterismIndices = self.currentFieldAsterismsIndices[astIndex]
        super().setAsterismData()
        self.firstConfigCall = False

            
    def translateTriangle(self, triangle, radius, angle):
        tt =  triangle.copy()
        tt[0,:] += radius*np.cos(angle/180*np.pi)
        tt[1,:] += radius*np.sin(angle/180*np.pi)
        return tt


    def generateTriangles(self, xrange, yrange, samples, photons=1900.0):
        self.nfields = 4
        self.nfieldsSizes = [10,10,10,100]
        for s in self.nfieldsSizes:
            self.cumAstSizes.append(s+self.cumAstSizes[-1])
        self.trianglesList = []
        self.polarTrianglesList = []
        maxRadius = np.sqrt(xrange**2+yrange**2)
        # the base shape is an equilateral trangle of a given radius
        standardRadius = maxRadius/6.0        
        standardAngles = [30, 150, 270]        
        standardRadi = [standardRadius, standardRadius, standardRadius]
        # triangle is a 3x2 np.array, containing the x,y coordinates of its veritces
        standardTriangle = np.vstack( (polarToCartesian(np.array( [standardRadi, standardAngles])), np.ones(3)*photons) )
        # the base position is this one
        standardPosRadius = maxRadius/12.0
        standardPosAngle = 45
        # A. dependency on position
        # A.1 dependency on radius (zenith)
        for rt in np.linspace(0, maxRadius*0.9, samples):
            triangle = self.translateTriangle(standardTriangle, rt, standardPosAngle)
            self.trianglesList.append(triangle)
            pTriangle = np.vstack( (cartesianToPolar(triangle[0:2,:]), np.ones(3)*photons) )                
            self.polarTrianglesList.append(pTriangle)        
        # A.2 dependency on angle (azimuth)
        for at in np.linspace(0, 180, samples):
            triangle = self.translateTriangle(standardTriangle, standardPosRadius, at)
            self.trianglesList.append(triangle)
            pTriangle = np.vstack( (cartesianToPolar(triangle[0:2,:]), np.ones(3)*photons) )                
            self.polarTrianglesList.append(pTriangle)
        # B. dependency on size
        for scale in np.linspace(0.5*standardRadius/samples, maxRadius*0.7, samples):
            triangle = standardTriangle.copy()
            ll = np.sqrt(standardTriangle[0,:]**2+standardTriangle[1,:]**2)
            triangle[0,:] = scale * standardTriangle[0,:]/ll
            triangle[1,:] = scale * standardTriangle[1,:]/ll
            self.trianglesList.append(triangle)
            pTriangle = np.vstack( (cartesianToPolar(triangle[0:2,:]), np.ones(3)*photons) )
            self.polarTrianglesList.append(pTriangle)
        # B. dependency on shape
        for a1 in np.linspace(1, 360, samples):
            for a2 in np.linspace(1, 180, samples):
                angles = [0, a1, a2]
                radii = standardRadi
                triangleS = np.vstack( (polarToCartesian(np.array( [radii, angles])), np.ones(3)*photons) )
                triangle = self.translateTriangle(triangleS, 3*standardPosRadius, standardPosRadius)
                self.trianglesList.append(triangle)                
                pTriangle = np.vstack( (cartesianToPolar(triangle[0:2,:]), np.ones(3)*photons) )                
                self.polarTrianglesList.append(pTriangle)
        return np.array(self.trianglesList), np.array(self.polarTrianglesList)


    def freqFromMagnitudeMAVIS(self, magnitude):
        return 900.0/(1 + np.exp((magnitude-17.5)*5.0)) + 100.0
#        if magnitude>19:
#            return 100.0
#        elif magnitude>18:
#            return 200.0
#        elif magnitude>17:
#            return 400.0
#        else:
#            return 1000.0


    def freqFromMagnitudeERIS(self, magnitude):
        if magnitude>17:
            return 100.0
        elif magnitude>16:
            return 250.0
        else:
            return 500.0


    def freqFromMagnitude(self, m):
        if 'freqRule' in self.my_data_map['ASTERISM_SELECTION'].keys():
            if self.my_data_map['ASTERISM_SELECTION']['freqRule'] == 'MORFEO':
                return 400.0/(1+np.exp((m-18.0)*2.5)) + 100.0
        else:
            if self.isMono:
                return self.freqFromMagnitudeERIS(m)
            else:
                return self.freqFromMagnitudeMAVIS(m)


    def freqsFromMagnitudes(self, magnitudes):
        result = []
        for m in magnitudes:
            result.append(self.freqFromMagnitude(m))
        return result


    def fluxFromMagnitude(self, m, b):
        m0 = 0.0
        f0 = None
        if 'flux'+b+'0' in self.my_data_map['ASTERISM_SELECTION'].keys():
            f0 = self.my_data_map['ASTERISM_SELECTION']['flux'+b+'0']
        else:
            if self.isMono:
                f0 = 1.51e10
            else:
                if b=='H':
                    f0 = 2.68e9
                if b=='J':
                    f0 = 3.72e9
        return f0 * np.power(10.0, (-(m-m0)/2.5) )


    def sampleMagnitudesAndFluxes(self):
        magnitudes = []
        fluxes = []
        for b in self.bands:
            if b=='H':
                mm = np.random.normal( (self.magnitudesRange[0]+self.magnitudesRange[1])/2.0, 3.0)
            elif b=='J':
                mm = magnitudes[-1] + np.random.uniform(0.5, 0.7)
            elif b=='R':
                mm = np.random.normal( (self.magnitudesRange[0]+self.magnitudesRange[1])/2.0, 3.0)
            elif b=='I':
                mm = magnitudes[-1] + np.random.uniform(0.5, 0.7)
            else:
                print('ERROR: Unexpected Band')
                mm = None
            magnitudes.append(mm)
            fluxes.append(self.fluxFromMagnitude(mm, b) ) # * self.freqFromMagnitude(mm)
        return magnitudes, fluxes

    def generateRandom(self, max_field):
        self.resetFieldsData()
        total_skipped_fields = 0
        total_skipped_asterisms = 0
        for i in range(self.globalOffset, self.globalOffset+max_field, 1):
            setsList = []
            polarSetsList = []
            # number_of_stars == number of asterisms when isMono
            number_of_stars0 = np.random.randint(self.minStars,self.maxStars)
            number_of_stars = number_of_stars0
            number_of_stars = 0
            self.reset_currentFieldsSourcesData()
            for j in range(number_of_stars0):
                xcoords = np.random.uniform(-self.fovRange[0],self.fovRange[0])
                ycoords = np.random.uniform(-self.fovRange[1],self.fovRange[1])
                flux = 0.0
                freq = 0
                magnitudes, fluxes = self.sampleMagnitudesAndFluxes()
                for b, mm, ff in zip(self.bands, magnitudes, fluxes):
                    freq = self.freqFromMagnitude(mm)
                    flux += ff * self.fluxScaling / freq
                if np.min(flux)<=0.0:
                    total_skipped_asterisms += 1
                    if not self.isMono:
                        number_of_stars -= 1
                    continue
                pcoords = cartesianToPolar( np.asarray([xcoords, ycoords]))
                source = np.array([pcoords[0], pcoords[1], flux, freq])
                self.appendSource(source)
                if self.isMono:
                    cartstar = np.vstack( [[xcoords], [ycoords], [flux], [freq]] )
                    polarstar = np.vstack( [ [pcoords[0]], [pcoords[1]], [flux], [freq] ])
                else:
                    cartstar = np.vstack( [xcoords, ycoords, flux, freq] )
                    polarstar = np.vstack( [pcoords[0], pcoords[1], flux, freq] )
                setsList.append(cartstar)
                polarSetsList.append(polarstar)
                number_of_stars += 1
            number_of_stars = max(0, number_of_stars)
            print('number_of_stars', number_of_stars)
            if self.isMono:
                number_of_asterisms = number_of_stars
                all_combos = list(itertools.combinations(list(range(number_of_stars)), 1))
            else:
                all_combos = list(itertools.combinations(list(range(number_of_stars)), 3))
                number_of_asterisms = len(all_combos)
            self.addFieldDataCombos(all_combos, number_of_asterisms, number_of_stars)
        print('\ntotal_skipped_fields: ', total_skipped_fields)
        print('total_skipped_asterisms: ', total_skipped_asterisms)
        print('total good asterisms: ', self.cumAstSizes[-1])
        print('total good asterisms: ', self.cumAstSizes)
        mString = ''
        if not self.isMono:
            mString = 'Multi'
        self.asterismsInputDataCartesian = np.array(self.setsList)
        self.asterismsInputDataPolar = np.array(self.polarSetsList)
        self.allAsterismsIndices = np.array(self.allAsterismsIndices)
        np.save(os.path.join(self.outputDir, self.file_field_simul + 'C.npy'), np.array( self.setsList))
        np.save(os.path.join(self.outputDir, self.file_field_simul + 'P.npy'), np.array( self.polarSetsList))
        np.save(os.path.join(self.outputDir, self.file_field_simul + 'F.npy'), np.array(self.nfieldsSizes))
        np.save(os.path.join(self.outputDir, self.file_field_simul + 'S.npy'), np.array(self.cumAstSizes))
        np.save(os.path.join(self.outputDir, self.file_field_simul + 'ST.npy'), np.array(self.cumStarSizes))
        np.save(os.path.join(self.outputDir, self.file_field_simul + 'IDX.npy'), self.allAsterismsIndices)


    def asterismDataFromRecArray(self, i, j):
        xcoords = self.asterismsRecArray[i][j]['COORD'][0]
        ycoords = self.asterismsRecArray[i][j]['COORD'][1]
        fluxes = np.zeros(len(xcoords))
        for b in self.bands:
            mm = self.asterismsRecArray[i][j][b+'MAG']
            freqs = self.freqsFromMagnitudes(mm)
            if 'flux'+b+'0' in self.my_data_map['ASTERISM_SELECTION'].keys():
                if np.min(mm) == 0:
                    mm[np.where(mm == 0)] = 30
                fluxes += np.asarray(self.fluxFromMagnitude(mm, b)* self.fluxScaling / np.asarray(freqs)) + 1e-3
            else:
                ff = self.asterismsRecArray[i][j]['FLUX' + b]
                fluxes += np.abs(np.asarray(ff) * self.fluxScaling / np.asarray(freqs)) + 1e-3
        return xcoords, ycoords, fluxes, freqs


    def generateFromRecArrayMulti(self, max_field):
        self.resetFieldsData()
        total_skipped_fields = 0
        total_skipped_asterisms = 0
        for i in range(self.globalOffset, self.globalOffset+max_field, 1):
            if self.progressStatus:
                sys.stdout.write('\r')
                sys.stdout.write('generateFromRecArrayMulti: ' + str(i+1) + ' of ' + str(max_field))
                sys.stdout.flush()
                time.sleep(0.001)
            setsList = []
            polarSetsList = []
            if self.verbose:
                print('Loading Field')
            skipped_field = False
            if type(self.asterismsRecArray[i]) is np.int16 or type(self.asterismsRecArray[i]) is np.int64:
                if self.verbose:
                    print("Field:" + str(i) + " SKIPPED")
                total_skipped_fields += 1
                self.nfields -=1
                skipped_field = True
                self.skippedFieldIndexes.append(i)
            else:
                if skipped_field:
                    number_of_asterisms0 = 0
                else:
                    number_of_asterisms0 = len(self.asterismsRecArray[i])
                    if self.verbose:
                        print('Potential number of asterisms', number_of_asterisms0)
                    number_of_asterisms = number_of_asterisms0
                    number_of_asterisms = 0
                    self.reset_currentFieldsSourcesData()
                    number_of_stars = 0
                    all_combos = []
                    for j in range(number_of_asterisms0):
                        xcoords, ycoords, fluxes, freqs = self.asterismDataFromRecArray(i, j)
                        #if np.min(fluxes)<=0.0:
                        #    print("Field:" + str(i) + "- Asterism:" + str(j) + " SKIPPED because of Flux 0 star")
                        #    total_skipped_asterisms += 1
                        #    number_of_asterisms -= 1
                        #    continue
                        s_index = []
                        for si in range(3):
                            pcoords = cartesianToPolar(np.asarray([xcoords[si], ycoords[si]]))
                            source = np.array([pcoords[0], pcoords[1], fluxes[si], freqs[si]])
                            ss = self.sourceIsPresent(source)
                            if ss==-1:
                                self.appendSource(source)
                                ss = number_of_stars
                                number_of_stars += 1
                            s_index.append(ss)
                        all_combos.extend([s_index])
                        pcoords = cartesianToPolar( np.asarray([xcoords, ycoords]))
                        asterism = np.vstack([xcoords, ycoords, fluxes, freqs])
                        pasterism = np.vstack([pcoords[0,:], pcoords[1,:], fluxes, freqs])
                        setsList.append(asterism)
                        polarSetsList.append(pasterism)
                number_of_asterisms = len(all_combos)
                number_of_asterisms = max(0, number_of_asterisms)
                if self.verbose:
                    print('number of asterisms', number_of_asterisms)
                if number_of_asterisms==0:
                    total_skipped_fields +=1
                    # self.cumAstSizes.append(self.cumAstSizes[-1])
                    # self.cumStarSizes.append(self.cumStarSizes[-1])
                    continue
                self.addFieldDataCombos(all_combos, number_of_asterisms, number_of_stars)
        print('\ntotal_skipped_fields: ', total_skipped_fields)
        print('total_skipped_asterisms: ', total_skipped_asterisms)
        print('total good asterisms: ', self.cumAstSizes[-1])
        self.asterismsInputDataCartesian = np.array(self.setsList)
        self.asterismsInputDataPolar = np.array(self.polarSetsList)
        self.allAsterismsIndices = np.array(self.allAsterismsIndices)


    def generateFromRecArray(self, max_field):
        self.resetFieldsData()
        total_skipped_fields = 0
        total_skipped_asterisms = 0
        for i in range(self.globalOffset, self.globalOffset+max_field, 1):
            if self.progressStatus:
                sys.stdout.write('\r')
                sys.stdout.write('generateFromRecArray: ' + str(i+1) + ' of ' + str(max_field))
                sys.stdout.flush()
                time.sleep(0.001)
            skipped_field = False
            if type(self.asterismsRecArray[i]) is np.int16 or type(self.asterismsRecArray[i]) is np.int64:
                if self.verbose:
                    print("Field:" + str(i) + " SKIPPED")
                total_skipped_fields += 1
                self.nfields -=1
                skipped_field = True
                self.skippedFieldIndexes.append(i)
            else:
                if skipped_field:
                    number_of_asterisms0 = 0
                else:
                    number_of_asterisms0 = len(self.asterismsRecArray[i])
                    if self.verbose:
                        print('Potential number of asterisms', number_of_asterisms0)
                    number_of_asterisms = number_of_asterisms0
                    number_of_asterisms = 0
                    number_of_stars = 0
                    self.reset_currentFieldsSourcesData()
                    for j in range(number_of_asterisms0):
                        xcoords, ycoords, fluxes, freqs = self.asterismDataFromRecArray(i, j)
                        if np.min(fluxes)<=0.0:
                            if self.verbose:
                                print("Field:" + str(i) + "- Asterism:" + str(j) + " SKIPPED because of Flux 0 star")
                            total_skipped_asterisms += 1                            
                            #if not self.isMono:
                            number_of_asterisms -= 1                            
                            continue
                        for si in range(3):
                            ccoords = np.asarray([xcoords[si], ycoords[si]] )
                            pcoords = cartesianToPolar( np.asarray([xcoords[si], ycoords[si]] ) )
                            source = np.array([pcoords[0], pcoords[1], fluxes[si], freqs[si]])
                            s_index = self.sourceIsPresent(source)
                            if s_index==-1 and np.abs(ccoords[0])<0.5*self.my_data_map['telescope']['TechnicalFoV'] and \
                               np.abs(ccoords[1])<0.5*self.my_data_map['telescope']['TechnicalFoV']:
                                self.appendSource(source)
                                asterism = np.vstack( [[xcoords[si]], [ycoords[si]], [fluxes[si]], [freqs[si]]] )
                                pasterism = np.vstack( [ [pcoords[0]], [pcoords[1]], [fluxes[si]], [freqs[si]]] )
                                self.setsList.append(asterism)
                                self.polarSetsList.append(pasterism)
                                number_of_asterisms += 1
                                number_of_stars += 1
                    number_of_asterisms = max(0, number_of_asterisms)
                    all_combos = list(itertools.combinations(list(range(number_of_asterisms)), 1))
                    self.updateAsterismIndices(all_combos, number_of_asterisms, number_of_stars)
        #    print('number_of_asterisms', number_of_asterisms)
        print('\ntotal_skipped_fields: ', total_skipped_fields)
        print('total_skipped_asterisms: ', total_skipped_asterisms)
        print('total good asterisms: ', self.cumAstSizes[-1])
        self.asterismsInputDataCartesian = np.array(self.setsList)
        self.asterismsInputDataPolar = np.array(self.polarSetsList)
        self.allAsterismsIndices = np.array(self.allAsterismsIndices)


    def sourceIsPresent(self, source):
        # a source is a triple (array) of zenith, azimuth, photons
        for z, a, n, f, ii in zip(self.currentFieldsSourcesData['Zenith'],
                               self.currentFieldsSourcesData['Azimuth'],
                               self.currentFieldsSourcesData['NumberPhotons'],
                               self.currentFieldsSourcesData['Frequencies'],
                               range(len(self.currentFieldsSourcesData['Zenith'])) ):
            if np.allclose(np.array([z,a,n,f]), source):
                return ii
        return -1


    def getSourcesData(self, fields):
        self.reset_currentFieldsSourcesData()
        for field in fields:
            fieldsize = self.nfieldsSizes[field]
            fieldsizeStars = self.cumStarSizes[field+1] - self.cumStarSizes[field]
            if self.verbose:
                print('fieldsize:', fieldsizeStars, ' Stars')
            firstStarInAsterismIndex = self.cumStarSizes[field+1]
            doneStars = {}
            for s in range(fieldsizeStars):
                doneStars[int(s)] = False
            for ast in range(fieldsize):
                if self.allAsterismsIndices.size==0:
                    continue
                astIndexGlobal = self.cumAstSizes[field]+ast
                astIndices = self.allAsterismsIndices[astIndexGlobal]
                self.currentFieldAsterismsIndices.append(astIndices)
                for s, si in zip(astIndices, range(3)):
                    if not doneStars[s]:
                        z = self.asterismsInputDataPolar[astIndexGlobal, 0, si]
                        a = self.asterismsInputDataPolar[astIndexGlobal, 1, si]
                        n = self.asterismsInputDataPolar[astIndexGlobal, 2, si]
                        f = self.asterismsInputDataPolar[astIndexGlobal, 3, si]
                        source = np.array([z, a, n, f])
                        self.appendSource(source)
                        doneStars[s] = True


    def selectData(self, fieldIndex1, fieldIndex2 = None):
        if fieldIndex2 is None:
            self.currentFieldsize = self.nfieldsSizes[fieldIndex1]
            self.currentField = fieldIndex1
            self.getSourcesData([fieldIndex1])
            self.currentBase = self.cumAstSizes[fieldIndex1]
            self.covsarray = np.array(self.cov_ellipses_Asterism)[self.currentBase:self.currentBase+self.currentFieldsize, 0,1]**2 + np.array(self.cov_ellipses_Asterism)[self.currentBase:self.currentBase+self.currentFieldsize, 0,2]**2
            self.strehls = np.array(self.strehl_Asterism)[self.currentBase:self.currentBase+self.currentFieldsize][:, 0]
            self.penalties = np.array(self.penalty_Asterism)[self.currentBase:self.currentBase+self.currentFieldsize][:, 0]
        else:
            self.currentFieldsize = 0
            fieldIndices = list(range(fieldIndex1, fieldIndex2, 1))
            for ii in fieldIndices:
                self.currentFieldsize += self.nfieldsSizes[ii]
            self.getSourcesData(fieldIndices)
            self.currentBase = self.cumAstSizes[fieldIndex1]
            self.covsarray = np.array(self.cov_ellipses_Asterism)[self.currentBase:self.currentBase+self.currentFieldsize, 0,1]**2 + np.array(self.cov_ellipses_Asterism)[self.currentBase:self.currentBase+self.currentFieldsize, 0,2]**2
            self.strehls = np.array(self.strehl_Asterism)[self.currentBase:self.currentBase+self.currentFieldsize][:, 0]
            self.penalties = np.array(self.penalty_Asterism)[self.currentBase:self.currentBase+self.currentFieldsize][:, 0]
        if self.covsarray.shape[0]!=0:
            self.jitterMeasure = self.penalties
            self.minJitter_id = np.argmin(self.jitterMeasure)+self.currentBase
            self.sortedJitterIndices = np.argsort(self.jitterMeasure, axis=0)
            self.lastJitterIndex = self.jitterMeasure.shape[0] + self.currentBase


    def fitHeuristicModel(self, fieldIndex1, fieldIndex2, modelName, num_epochs=500, steps = [1e-4], geom=[128]*8):
        self.selectData(fieldIndex1, fieldIndex2)
        self.fitModel(modelName, num_epochs, steps, geom)


    def setModelData(self):
        self.jitterM = np.log(self.jitterMeasure+1)
        self.coordsM = self.asterismsInputDataCartesian
        self.xcoordsM = self.coordsM[self.currentBase:self.lastJitterIndex, 0, :]
        self.ycoordsM = self.coordsM[self.currentBase:self.lastJitterIndex, 1, :]
        self.rcoordsM = self.asterismsInputDataPolar[self.currentBase:self.lastJitterIndex, 0, :]
        self.freqsM  = self.asterismsInputDataCartesian[self.currentBase:self.lastJitterIndex, 3, :]
        self.fluxesM = self.asterismsInputDataCartesian[self.currentBase:self.lastJitterIndex, 2, :]
        self.fluxesM  = np.log(self.fluxesM)
        if not self.isMono:
            self.inputDataT = torch.hstack((torch.from_numpy(np.sqrt(self.xcoordsM**2+self.ycoordsM**2)).double().clone().detach().requires_grad_(True),
                               torch.from_numpy(np.arctan2(self.xcoordsM, self.ycoordsM)).double().clone().detach().requires_grad_(True),
                               torch.from_numpy(self.fluxesM).double().clone().detach().requires_grad_(True),
                               torch.from_numpy(self.freqsM).double().clone().detach().requires_grad_(True)) )
            self.jitterT = torch.from_numpy(self.jitterM).double().clone().detach().requires_grad_(True)
            self.jitterT = self.jitterT[:, None].requires_grad_(True)


    def fitModel(self, modelName, num_epochs, steps, geom):
        if self.isMono:
            ast = 0
            self.firstConfigCall = True
            self.currentField = 0
            self.getSourcesData([0])
            self.currentAsterism = ast
            self.doOverallSimulation(ast)
        self.setModelData()
        if self.isMono:
            trainInput = np.abs(np.array( [self.rcoordsM[:,0], self.fluxesM[:,0]] )) # /self.freqsM[:,0]
            # 100, 250, 500
            idxF0 = np.where(self.freqsM[:,0]<200)
            idxF1 = np.where((self.freqsM[:,0]<350) & (self.freqsM[:,0]>=200))
            idxF2 = np.where(self.freqsM[:,0]>=350)
            idxV = [idxF0, idxF1, idxF2]
            jitterTrainM = np.abs(self.jitterM)
            jitterTrain = np.exp(jitterTrainM)-1
            # func = funcPolar
            # popt, pcov = curve_fit(func, trainInput, jitterTrain)
            self.monoModel = []
            for idx in idxV:
                ww = np.power(1/(jitterTrain[idx]), 2)
                smoothing = 10*jitterTrain[idx].shape[0]
                self.monoModel.append(interpolate.SmoothBivariateSpline(trainInput[0, idx], trainInput[1,idx], jitterTrainM[idx], w=ww, kx=4, ky=4))
#            monoModel = interpolate.SmoothBivariateSpline(trainInput[0,:], trainInput[1,:], jitterTrainM, w=ww, kx=5, ky=5)
            self.monoModel.append(self.pointings_FWHM_mas)
            self.monoModel.append(self.HO_res)
            # print(popt, pcov)
            # for i, j in zip(popt, ascii_uppercase):
            #     print(f"{j} = {i:.6f}")
            #jitterApproxTrain = monoModel(trainInput, *popt)
            jitterApproxTrainM = np.zeros(trainInput.shape[1])            
            for i, idx in enumerate(idxV):
                grid_x = trainInput[0,idx]
                grid_y = trainInput[1,idx]
                jitterApproxTrainM[idx] = self.monoModel[i].__call__(grid_x, grid_y, grid=False)
            current_pointings_FWHM_mas = self.monoModel[len(idxV)]
            current_HO_res = np.asarray(self.monoModel[len(idxV)+1].get())
            
            jitterApproxTrain = np.exp(jitterApproxTrainM)-1

            print(jitterApproxTrain.shape, jitterTrain.shape)
            absoluteErrorTrain = np.abs((jitterTrain-jitterApproxTrain))
            print( "Mean Absolute Error Train", np.mean(absoluteErrorTrain))
            relativeErrorTrain = 2.0 * np.abs((jitterTrain-jitterApproxTrain)/np.abs(jitterTrain+jitterApproxTrain))
            print( "Mean Relative Error Train", np.mean(relativeErrorTrain))
            absoluteErrorTrain = np.abs((jitterTrain-jitterApproxTrain))
            print( "Median Absolute Error Train", np.median(absoluteErrorTrain))
            rmsErrorTrain = np.sqrt(np.mean( (jitterTrain-jitterApproxTrain)*(jitterTrain-jitterApproxTrain) ) )
            print( "RMS Error Train", rmsErrorTrain)
            if self.doPlotAst:
                plt.figure(figsize=(10, 10))
                plt.scatter(trainInput[0,:], absoluteErrorTrain, alpha=0.2, s=4)
                plt.xlabel('NGS distance')
                plt.ylabel('absolute error')
                plt.yscale('log')
                plt.show()
                plt.figure(figsize=(10, 10))
                plt.scatter(trainInput[1,:], absoluteErrorTrain, alpha=0.2, s=4)
                plt.xlabel('flux')
                plt.ylabel('absolute error')
                plt.xscale('log')
                plt.yscale('log')
                plt.show()
                plt.figure(figsize=(10, 10))
                plt.scatter(jitterTrain, absoluteErrorTrain, alpha=0.2, s=4)
                plt.xlabel('penalty')
                plt.ylabel('absolute error')
                plt.yscale('log')
                plt.show()
        else:
            modelNameFile = os.path.join(self.outputDir, modelName + '.pth')
            model = NeuralNetwork(12, geom)
            model.setData(self.inputDataT, self.jitterT, 0.1, True)
            trainInput, jitterTrainM, jitterApproxTrainM, testInput, jitterTestM, jitterApproxTestM = model.trainModel(num_epochs, steps, modelNameFile)

            jitterTrain = np.exp(jitterTrainM)-1
            jitterApproxTrain = np.exp(jitterApproxTrainM)-1
            jitterTest = np.exp(jitterTestM)-1
            jitterApproxTest = np.exp(jitterApproxTestM)-1

            absoluteErrorTrain = np.abs((jitterTrain-jitterApproxTrain))
            print( "Mean Absolute Error Train", np.mean(absoluteErrorTrain))
            absoluteErrorTest = np.abs((jitterTest-jitterApproxTest))
            print( "Mean Absolute Error Test", np.mean(absoluteErrorTest))
            absoluteErrorTrain = np.abs((jitterTrain-jitterApproxTrain))
            print( "Median Absolute Error Train", np.median(absoluteErrorTrain))
            absoluteErrorTest = np.abs((jitterTest-jitterApproxTest))
            print( "Median Absolute Error Test", np.median(absoluteErrorTest))
            rmsErrorTrain = np.sqrt(np.mean( (jitterTrain-jitterApproxTrain)*(jitterTrain-jitterApproxTrain) ) )
            print( "RMS Error Train", rmsErrorTrain)
            rmsErrorTest = np.sqrt(np.mean( (jitterTest-jitterApproxTest)*(jitterTest-jitterApproxTest) ) )
            print( "RMS Error Test", absoluteErrorTest)
            if self.doPlotAst:
                plt.figure(figsize=(10, 10))
                plt.scatter(trainInput[:,0], absoluteErrorTrain, alpha=0.2, s=4)
                plt.title('input 0, train set')
                plt.ylabel('absolute error')
                plt.show()
                plt.figure(figsize=(10, 10))
                plt.scatter(trainInput[:,3], absoluteErrorTrain, alpha=0.2, s=4)
                plt.title('input 3, train set')
                plt.ylabel('absolute error')
                plt.show()
                plt.figure(figsize=(10, 10))
                plt.scatter(trainInput[:,6], absoluteErrorTrain, alpha=0.2, s=4)
                plt.title('input 6, train set')
                plt.ylabel('absolute error')
                plt.show()
                plt.figure(figsize=(10, 10))
                plt.scatter(trainInput[:,9], absoluteErrorTrain, alpha=0.2, s=4)
                plt.title('input 9, train set')
                plt.ylabel('absolute error')
                plt.show()
                plt.figure(figsize=(10, 10))
                plt.scatter(jitterTrain, absoluteErrorTrain, alpha=0.2, s=4)
                plt.title('penalty, train set')
                plt.ylabel('absolute error')
                plt.show()
                plt.figure(figsize=(10, 10))
                plt.scatter(testInput[:,0], absoluteErrorTest, alpha=0.2, s=4)
                plt.title('input 0, test set')
                plt.ylabel('absolute error')
                plt.show()
                plt.figure(figsize=(10, 10))
                plt.scatter(testInput[:,3], absoluteErrorTest, alpha=0.2, s=4)
                plt.title('input 3, test set')
                plt.ylabel('absolute error')
                plt.show()
                plt.figure(figsize=(10, 10))
                plt.scatter(testInput[:,6], absoluteErrorTest, alpha=0.2, s=4)
                plt.title('input 6, test set')
                plt.ylabel('absolute error')
                plt.show()
                plt.figure(figsize=(10, 10))
                plt.scatter(testInput[:,9], absoluteErrorTest, alpha=0.2, s=4)
                plt.title('input 9, test set')
                plt.ylabel('absolute error')
                plt.show()
                plt.figure(figsize=(10, 10))
                plt.scatter(jitterTest, absoluteErrorTest, alpha=0.2, s=4)
                plt.title('penalty, test set')
                plt.ylabel('absolute error')
                plt.show()
        if self.isMono:
#            np.save(os.path.join(self.outputDir, modelName +'.npy'), np.array(popt))
            with open(os.path.join(self.outputDir, modelName)+'.npy', 'wb') as f:
                pickle.dump(self.monoModel, f)
#        else:
#            model.save_model(os.path.join(self.outputDir, modelName + '.pth'))


    def runHeuristicModel(self):
        # if self.isMono: # only usable for mono for now
        self.rcoordsM = self.asterismsInputDataPolar[:, 0, :]
        self.fluxesM = self.asterismsInputDataPolar[:, 2, :]        
        self.freqsM = self.asterismsInputDataPolar[:, 3, :]
#        monoModel = funcPolar
        inputDataTestCpu = np.abs(np.array( [self.rcoordsM[:,0], self.fluxesM[:,0]] )) # /self.freqsM[:,0]
#        jitterApprox = monoModel(inputDataTestCpu, *self.monoModel)
        # 100, 250, 500
        idxF0 = np.where(self.freqsM[:,0]<200)
        idxF1 = np.where((self.freqsM[:,0]<350) & (self.freqsM[:,0]>=200))
        idxF2 = np.where(self.freqsM[:,0]>=350)
        idxV = [idxF0, idxF1, idxF2]
        jitterApproxM = np.zeros(inputDataTestCpu.shape[1])            
        for i, idx in enumerate(idxV):
            if len(idx)>0:
                grid_x = inputDataTestCpu[0,idx]
                grid_y = inputDataTestCpu[1,idx]
                jitterApproxM[idx] = self.monoModel[i].__call__(grid_x, grid_y, grid=False)
        current_pointings_FWHM_mas = np.asarray(self.monoModel[len(idxV)])
        current_HO_res = np.asarray(self.monoModel[len(idxV)+1].get())
        jitterApprox = np.exp(jitterApproxM)-1
        strehls = np.exp( -4*np.pi**2 * ( (jitterApprox)**2 )/(self.wvl*1e9)**2)
        vv = jitterApprox**2 - current_HO_res**2
        vv[np.where(vv<0)] = 0
        lo_res = np.sqrt(vv)
        scale = (np.pi/(180*3600*1000) * self.TelescopeDiameter / (4*1e-9))
        fwhms_lo = 2.355 * lo_res/scale / np.sqrt(2)
        fwhms = np.sqrt(fwhms_lo**2 + current_pointings_FWHM_mas**2)
        ees = [0]*fwhms.shape[0]
        astList = self.asterismsInputDataPolar.tolist()
        sortedJitterIndicesModel = np.argsort(jitterApprox, axis=0)
        jitterApproxSorted = jitterApprox[sortedJitterIndicesModel]
        deltas = np.abs( (jitterApproxSorted[:-1] - jitterApproxSorted[1:])/jitterApproxSorted[:-1] )
        cumDeltas = np.cumsum(deltas)
        ss = 1
        Un = 0.12
        if deltas[0]<Un:
            ss = np.where(cumDeltas<Un)[0].shape[0]+1
            print('Low Delta between best and n-th best asterism, n:', ss, 'of ', jitterApprox.shape[0])
        return self.assembleOtuput(astList, sortedJitterIndicesModel.tolist(), jitterApprox.tolist(), strehls.tolist(), fwhms.tolist(), ees)

    def assembleOtuput(self, astList, sortedJitterIndicesModel, jitterApprox, strehls, fwhms, ees, absolute=False):
        results = []
        for ii, index in enumerate(sortedJitterIndicesModel):
            stars = astList[index]
            zeniths = stars[0]
            azimuths = stars[1]
            photons = stars[2]
            freqs = stars[3]
            asterism = []
            for z, a, p, f in zip(zeniths, azimuths, photons, freqs):
                asterism.append(Star(z,a,p,f))
            if absolute:
                resultElement = AsterismProperties(index, asterism, jitterApprox[ii], strehls[ii], fwhms[ii], ees[ii])
            else:
                resultElement = AsterismProperties(index, asterism, jitterApprox[index], strehls[index], fwhms[index], ees[index])
            results.append(resultElement)
        return results

    def testHeuristicModel(self, fieldIndex1, fieldIndex2, modelName, geom):
        self.selectData(fieldIndex1, fieldIndex2)
        self.setModelData()
        func = None
        model = None
        if self.isMono:
            inputIndicesPlots = [0,1]
            with open(os.path.join(self.outputDir, modelName)+'.npy', 'rb') as f:
                popt = pickle.load(f)
#            popt = np.load(os.path.join(self.outputDir, modelName)+'.npy')                        
#            func = funcPolar
            funcbs = popt            
            inputDataTestCpu = np.abs(np.array( [self.rcoordsM[:,0], self.fluxesM[:,0]] ))            
            jitterTestCpuM = np.abs(self.jitterM)
            # jitterApprox = monoModel(inputDataTestCpu, *popt)
            # 100, 250, 500
            idxF0 = np.where(self.freqsM[:,0]<200)
            idxF1 = np.where((self.freqsM[:,0]<350) & (self.freqsM[:,0]>=200))
            idxF2 = np.where(self.freqsM[:,0]>=350)
            idxV = [idxF0, idxF1, idxF2]
            jitterApproxM = np.zeros(inputDataTestCpu.shape[1])
            for i, idx in enumerate(idxV):
                grid_x = inputDataTestCpu[0,idx]
                grid_y = inputDataTestCpu[1,idx]
                jitterApproxM[idx] = self.monoModel[i].__call__(grid_x, grid_y, grid=False)
            current_pointings_FWHM_mas = self.monoModel[len(idxV)]
            current_HO_res = np.asarray(self.monoModel[len(idxV)+1].get())
            inputDataTestCpu = inputDataTestCpu.transpose()
        else:
            inputIndicesPlots = [0,3,6,8]
            model_path = os.path.join(self.outputDir, modelName +'.pth')
            model = NeuralNetwork(12, geom)
            model.load_model(model_path)
            model.to(device)
            model.setData(self.inputDataT, self.jitterT, 1.0, True)
            model.eval()
            inputDataTest, jitterTestM = model.test_loader.dataset[:]
            with torch.no_grad():
                jitterApproxT = model(inputDataTest)
            jitterApproxM = jitterApproxT.detach().cpu().numpy()
            jitterApproxM = jitterApproxM[:,0]
            ## from here should be common 
            jitterTestCpuM =  jitterTestM[:,0].detach().cpu().numpy()
            inputDataTestCpu = inputDataTest.detach().cpu().numpy()

        jitterTestCpu = np.exp(jitterTestCpuM) - 1
        jitterApprox = np.exp(jitterApproxM) - 1

        signedError = (jitterTestCpu-jitterApprox)
        absoluteError = np.abs(signedError)
        relativeError = absoluteError/jitterTestCpu
        sortedJitterIndicesModel = np.argsort(jitterApprox, axis=0)
        rmsErrorTest = np.sqrt(np.mean( np.where(relativeError<1.0, relativeError*relativeError, 0)))
        Un = 1.6 * rmsErrorTest
        print( "Average Absolute Error", np.mean(absoluteError))
        print( "STD Absolute Error", np.std(absoluteError))
        print( "Average Relative Error", np.mean(relativeError))
        print( "STD Relative Error", np.std(relativeError))
        print( "RMS Error Test", rmsErrorTest)
        print("Un:", Un)
        # Un = np.mean(absoluteError) + 2 * np.std(absoluteError)
        if self.doPlotAst:
            for inputIndex in inputIndicesPlots:
                plt.figure(figsize=(10, 10))
                plt.scatter(inputDataTestCpu[:,inputIndex], relativeError, alpha=0.2, s=4)
                plt.ylabel('relative error')
                plt.xscale('log')
                plt.yscale('log')
                plt.show()
            plt.figure(figsize=(10, 10))
            plt.scatter(jitterTestCpu, relativeError, alpha=0.2, s=4)
            plt.xlabel('penalty')
            plt.ylabel('relative error')
            plt.yscale('log')
            plt.show()
        totalS = 0
        lv = signedError/jitterTestCpu
        sigma = np.std(lv)
        num_bins = 40
        if self.doPlotAst:
            fig = plt.figure(figsize=(10, 6))
            ax2 = fig.add_subplot(1, 1, 1)
            ax2.hist( lv.ravel(), bins=np.linspace(-3*sigma, 3*sigma, num=num_bins))
            plt.xlabel('relative error')
            plt.ylabel('counts')
            plt.show()
#        ax2.set_ylabel('Simulations')
#        ax2.set_xlabel('Error [' +  ts.data_loader.units[l_ind] + ']')
        totalAsterisms = 0
        bestIndexPlot = np.zeros(50)
        for field in range(fieldIndex1, fieldIndex2, 1):
            self.currentField = field
            self.getSourcesData([field])
            fieldsize = len(self.currentFieldAsterismsIndices)
            base = self.cumAstSizes[field]
            self.selectData(field)
            self.setModelData()
            if self.isMono:
                inputDataTest = np.abs(np.array( [self.rcoordsM[:,0], self.fluxesM[:,0]] ))            
                jitterTestCpu = np.abs(self.jitterM)            
#                jitterApprox = func(inputDataTest, *popt)
                # 100, 250, 500
                idxF0 = np.where(self.freqsM[:,0]<200)
                idxF1 = np.where((self.freqsM[:,0]<350) & (self.freqsM[:,0]>=200))
                idxF2 = np.where(self.freqsM[:,0]>=350)
                idxV = [idxF0, idxF1, idxF2]
                jitterApprox = np.zeros(inputDataTest.shape[1])            
                for i, idx in enumerate(idxV):
                    if len(idx)>0:
                        grid_x = inputDataTest[0,idx]
                        grid_y = inputDataTest[1,idx]
                        jitterApprox[idx] = self.monoModel[i].__call__(grid_x, grid_y, grid=False)
                inputDataTestCpu = inputDataTest.transpose()
            else:
                model.setData(self.inputDataT, self.jitterT, 1.0, False)
                inputDataTest, jitterTest = model.test_loader.dataset[:]
                with torch.no_grad():
                    jitterApproxT = model(inputDataTest)
                jitterApprox = jitterApproxT.detach().cpu().numpy()
                jitterApprox = jitterApprox[:,0]
                jitterTestCpu =  jitterTest[:,0].detach().cpu().numpy()
            if jitterApprox.shape[0]>1:
                absoluteError = np.abs((jitterTestCpu-jitterApprox)/jitterTestCpu)
                totalAsterisms += absoluteError.shape[0]
                sortedJitterIndicesModel = np.argsort(jitterApprox, axis=0)
                jitterApproxSorted = jitterApprox[sortedJitterIndicesModel]
                deltas = np.abs( (jitterApproxSorted[:-1] - jitterApproxSorted[1:])/jitterApproxSorted[:-1] )
                cumDeltas = np.cumsum(deltas)
                if deltas[0]<Un:
                    ss = np.where(cumDeltas<Un)[0].shape[0]+1
                    print(field, 'Low Delta between best and n-th best asterism, n:', ss, 'of ', absoluteError.shape[0])
                    totalS += ss
                else:
                    totalS += 1
                correctBestIndex = self.sortedJitterIndices[0]
                approxBestIndex = sortedJitterIndicesModel[0]
                bestIndexInApproxArray = np.where(sortedJitterIndicesModel == correctBestIndex)[0][0]
                bestIndexPlot[bestIndexInApproxArray] += 1.0
                if correctBestIndex != approxBestIndex:
                    print('Wrong order would be detected relying only on the model:')
                    print(field, deltas)
                    print("correct indices", self.sortedJitterIndices)
                    print("approx indices", sortedJitterIndicesModel)
                    print("approx jitter", jitterApproxSorted)
                    print("deltas", deltas)
                    print("cumDeltas", cumDeltas)
                    print('bestIndexInApproxArray', bestIndexInApproxArray)
                    if cumDeltas[bestIndexInApproxArray-1]>Un and bestIndexInApproxArray>2:
                        print('ERROR: Correct ordering estimation missed!!!!!!!!')
        ii = np.max(np.nonzero(bestIndexPlot))
        r = bestIndexPlot[:ii+2]
        r /= np.sum(r) / 100
        if self.doPlotAst:
            fig = plt.figure(figsize=(16, 8))
            ax = fig.add_subplot(1, 1, 1)
            ax.tick_params(axis='both', which='major', labelsize=20)
            idx = range(len(r))
            ax.bar( idx, r )
            ax.set_xticks(idx, labels=map(str, idx), fontsize=20)
            plt.xlabel('Postion of the estimated best asterism', fontsize=24)
            plt.ylabel('Times [%]', fontsize=24)
            ax.set_yscale('linear')
            ax.yaxis.grid(True)
            plt.show()
            print(bestIndexPlot)
            print('Number of asterisms having a rating too close to the best one', totalS, 'of ', totalAsterisms, 100 * totalS/totalAsterisms)

    def plotFieldInterval(self, fieldIndex1, fieldIndex2):
        self.selectData(fieldIndex1, fieldIndex2)
        print("*")
        print("* Plotting Fields " + str(fieldIndex1) + ' to ' + str(fieldIndex2) + " - number of asterisms :" + str(self.currentFieldsize))
        print("*")
        print('self.covsarray.shape', self.covsarray.shape)
        print('self.cov_ellipses_Asterism.shape', self.cov_ellipses_Asterism.shape)
        if self.covsarray.shape[0]!=0:
            self.twoPlots()


    def plotField(self, fieldIndex1):
        self.selectData(fieldIndex1)
        print("*")
        print("* Plotting Field " + str(fieldIndex1) + " - number of asterisms :" + str(self.currentFieldsize))
        print("*")
        print('self.covsarray.shape', self.covsarray.shape)
        print('len(self.cov_ellipses_Asterism)', len(self.cov_ellipses_Asterism))
        if self.covsarray.shape[0]!=0:
            self.twoPlots()


    def computeAsterisms(self, eeRadiusInMas, index=None, doConvolve=False, plotGS=False):
        if index==None:
            singleAsterism = False
        else:
            singleAsterism = True
        self.doConvolveAsterism = doConvolve
        self.eeRadiusInMas = eeRadiusInMas
        self.fwhm_Asterism = []
        self.ee_Asterism = []
        self.cov_ellipses_Asterism = []
        self.strehl_Asterism = []
        self.penalty_Asterism = []
        nf = self.nfields
        if singleAsterism:
            nf=1
        for field in range( min( nf, len(self.nfieldsSizes) )):
            if self.progressStatus:
                sys.stdout.write('\r')
                sys.stdout.write('computeAsterisms: ' + str(field) + ' of ' + str(min( nf, len(self.nfieldsSizes) )) )
                sys.stdout.flush()
                time.sleep(0.001)

            self.firstConfigCall = True
            if not singleAsterism:
                self.firstSimCall = True
            self.currentField = field
            if self.verbose:
                print('self.currentField:', self.currentField)
            if self.verbose:
                print('field', field)
            self.getSourcesData([field])
            fieldsize = len(self.currentFieldAsterismsIndices)
            if fieldsize==0:
                continue
            base = self.cumAstSizes[field]
            listOfAsterisms = [index]
            if index is None:
                 listOfAsterisms = list(range(fieldsize))
            if len(listOfAsterisms)==0:
                if self.verbose:
                    print('Skipping')
                continue
            if plotGS:
                self.plot_directions(base, len(listOfAsterisms))
            for ast in listOfAsterisms:
                self.currentAsterism = ast
#                try:
                self.doOverallSimulation(ast)
                self.computeMetrics()
#               except:
#                    print("Unexpected Error Computing (field, asterism): ", self.currentField, self.currentAsterism)
#                    print(self.currentFieldsSourcesData)
                self.strehl_Asterism.append(np.array( [cpuArray(x) for x in self.sr]))
                self.penalty_Asterism.append(np.array( [cpuArray(x) for x in self.penalty]))
                self.fwhm_Asterism.append(self.fwhm)
                self.ee_Asterism.append(self.ee)
                self.cov_ellipses_Asterism.append(self.cov_ellipses)
            if (field+1) % 10 == 0 and not singleAsterism:
                if self.verbose:
                    print("Field " + str(field) + " DONE")
                np.save(os.path.join(self.outputDir, self.simulName+'fw.npy'), np.array(self.fwhm_Asterism))
                np.save(os.path.join(self.outputDir, self.simulName+'ee.npy'), np.array(self.ee_Asterism))
                np.save(os.path.join(self.outputDir, self.simulName+'covs.npy'), np.array(self.cov_ellipses_Asterism))
                np.save(os.path.join(self.outputDir, self.simulName+'sr.npy'), np.array(self.strehl_Asterism))
                np.save(os.path.join(self.outputDir, self.simulName+'penalty.npy'), np.array(self.penalty_Asterism))
        if not singleAsterism:
            np.save(os.path.join(self.outputDir, self.simulName+'fw.npy'), np.array(self.fwhm_Asterism))
            np.save(os.path.join(self.outputDir, self.simulName+'ee.npy'), np.array(self.ee_Asterism))
            np.save(os.path.join(self.outputDir, self.simulName+'covs.npy'), np.array(self.cov_ellipses_Asterism))
            np.save(os.path.join(self.outputDir, self.simulName+'sr.npy'), np.array(self.strehl_Asterism))
            np.save(os.path.join(self.outputDir, self.simulName+'penalty.npy'), np.array(self.penalty_Asterism))

        if nf==1:
            self.firstConfigCall = True
            #    print('Actual penalt for asterism', ii, ':', np.log(simulation.penalty_Asterism[0][0]+1))
            #    print('Actual Strehel for asterism', ii, ':', simulation.strehl_Asterism)
            #    print('Actual FWHM for asterism', ii, ':', simulation.fwhm_Asterism)
            if singleAsterism:
                jitter = self.penalty_Asterism[:][:]
                astList = self.asterismsInputDataPolar.tolist()
                return self.assembleOtuput( astList, [index], jitter[0], self.strehl_Asterism[0],
                                            self.fwhm_Asterism[0], self.ee_Asterism[0], absolute=True )
            else:
                jitter = self.penalty_Asterism[:][:]
                sortedJitterIndices = np.argsort(jitter, axis=0).reshape(-1).tolist()
                astList = self.asterismsInputDataPolar.tolist()
                jitter = np.asarray(jitter).reshape(-1).tolist()
                strehl_Asterism = np.asarray(self.strehl_Asterism).reshape(-1).tolist()
                fwhm_Asterism = np.asarray(self.fwhm_Asterism).reshape(-1).tolist()
                ee_Asterism = np.asarray(self.ee_Asterism).reshape(-1).tolist()
                return self.assembleOtuput( astList, sortedJitterIndices, jitter, strehl_Asterism,
                                            fwhm_Asterism, ee_Asterism)

    def reloadResults(self):
        self.fwhm_Asterism = np.load(os.path.join(self.outputDir, self.simulName+'fw.npy'))
        self.ee_Asterism = np.load(os.path.join(self.outputDir, self.simulName+'ee.npy'))
        self.cov_ellipses_Asterism = np.load(os.path.join(self.outputDir, self.simulName+'covs.npy'))
        self.strehl_Asterism = np.load(os.path.join(self.outputDir, self.simulName+'sr.npy'))
        self.penalty_Asterism = np.load(os.path.join(self.outputDir, self.simulName+'penalty.npy'))


    def plotAsterisms(self):
        al = self.jitterMeasure
        al = (np.max(al) - al)/np.max(al)
        np.random.seed(12345)
        if not self.currentBase:
            self.currentBase = 0
        X = self.asterismsInputDataCartesian
        xcoords = X[self.currentBase:self.lastJitterIndex, 0, :]
        ycoords = X[self.currentBase:self.lastJitterIndex, 1, :]
        fluxes  = X[self.currentBase:self.lastJitterIndex, 2, :]
        max_flux = np.max(fluxes)
        scales = 0.5 * np.log(fluxes+np.exp(1.0))
        ntriangles = X.shape[0]
        fig = plt.figure(figsize=(10, 10), dpi=90)
        ax = fig.add_subplot(1,1,1)
        ax.axvline(0, c='black', linewidth=0.5)
        ax.axhline(0, c='black', linewidth=0.5)
        
        norm = matplotlib.colors.Normalize()
        norm.autoscale(al)
        cm1 = cm.get_cmap('rainbow')
        sm = matplotlib.cm.ScalarMappable(cmap=cm1, norm=norm)
        sm.set_array(al)
        
        for i in range(self.lastJitterIndex-1, self.currentBase-1, -1):
            jj = self.sortedJitterIndices[i-self.currentBase]
            aa = 1.0 # al[i-self.currentBase] * 0.9 + 0.1
            coords = np.transpose(X[jj+self.currentBase, :2, :])
            xx = np.sum(X[jj+self.currentBase, 0, :])/float(X.shape[2])
            yy = np.sum(X[jj+self.currentBase, 1, :])/float(X.shape[2])
            px = X[jj+self.currentBase, 0, :]
            py = X[jj+self.currentBase, 1, :]
            if X.shape[2]==1:
                circle1 = plt.Circle((xx, yy), scales[jj], color=cm1(norm(al[jj])), fill=False)
                ax.add_patch(circle1)
            else:
                ppp=ax.quiver([xx, xx, xx], [yy, yy, yy], px-xx, py-yy, color=cm1(norm(al[jj])), width=0.003, scale_units='xy', scale=1, alpha = aa )

        # draw best asterism
        coords = np.transpose(X[self.minJitter_id, :2, :])
        if X.shape[2]==1:
            circle1 = plt.Circle((coords[0,0], coords[0,1]), scales[self.minJitter_id-self.currentBase], color='r', fill=True, alpha=0.5)
            ax.add_patch(circle1)
        else:
            t1 = plt.Polygon(coords, alpha = 0.3, color='r')
            ax.add_patch(t1)

        ax.scatter(xcoords, ycoords, s=scales**2 * 50, c='yellow', edgecolors='y', marker='*', label='Natural Guide Stars')

        # draw science targets
        ax.scatter(self.xxSciencePointigs, self.yySciencePointigs, c='blue', marker='x', label='Science Targets')

        # draw LGSs
        zenithSrc  = self.my_data_map['sources_HO']['Zenith']
        azimuthSrc = self.my_data_map['sources_HO']['Azimuth']
        pointings = polarToCartesian(np.array( [zenithSrc, azimuthSrc]))
        self.xxLGSPointigs         = pointings[0,:]
        self.yyLGSPointigs         = pointings[1,:]
        ax.scatter(self.xxLGSPointigs, self.yyLGSPointigs, s= 100.0, c='green', marker='*', edgecolors='green', label='Laser Guide Stars')

        # Finalize plot
        ax.legend()
        ax.set_aspect('equal', adjustable='box')
        xylim = 0.5*self.my_data_map['telescope']['TechnicalFoV']
        ax.set_xlim([-xylim,xylim])
        ax.set_ylim([-xylim,xylim])
        cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
        cbar = fig.colorbar(sm, cax=cax) # Similar to fig.colorbar(im, cax = cax)
        cbar.set_label('Asterism Rating')
        plt.show()


    def twoPlots(self):
        plt.plot(self.jitterMeasure[self.sortedJitterIndices])
        plt.show()
        self.plotAsterisms()

    def plot_directions(self, base, nAsterisms):
        '''
        Polar plot with science and GS (sources_HO and sources_LO) directions
        '''

        ticks_interval = 5

        # SCIENCE
        th_sci = np.array(self.my_data_map['sources_science']['Azimuth'])
        rr_sci = np.array(self.my_data_map['sources_science']['Zenith'])
        th_sci = th_sci/180*np.pi
        fig = plt.figure('SOURCES DIRECTIONS', figsize=(5,5))
        ax = fig.add_subplot(111, polar=True)
        ax.tick_params(labelsize=10)
        ax.set_rlabel_position(225) # theta position of the radius labels
        ax.scatter(th_sci, rr_sci, marker='*', color='blue', s=120, label='sources_science')

        # LGS
        th_HO = np.array(self.my_data_map['sources_HO']['Azimuth'])
        rr_HO = np.array(self.my_data_map['sources_HO']['Zenith'])
        th_HO = th_HO/180*np.pi
        ax.scatter(th_HO, rr_HO, marker='*', color='green', s=120, label='sources_HO')

        # NGS
        rr_LO = self.asterismsInputDataPolar[base:base+nAsterisms, 0, :]
        th_LO = self.asterismsInputDataPolar[base:base+nAsterisms, 1, :]/180*np.pi
        fluxes = self.asterismsInputDataPolar[base:base+nAsterisms, 2, :]

        rr_LO = rr_LO.flatten()
        th_LO = th_LO.flatten()
        fluxes = fluxes.flatten()

        dummy, indices = np.unique(rr_LO, return_index=True)
        rr_LO = rr_LO[indices]
        th_LO = th_LO[indices]
        fluxes = fluxes[indices]

        max_flux = np.max(fluxes)
        scales = 0.5 * np.log(fluxes+np.exp(1.0))
        ax.scatter(th_LO, rr_LO, marker='*', color='red', s=120, label='sources_LO')

        # Set ticks position
        max_pos = np.max([rr_sci.max(), rr_HO.max(), rr_LO.max()]) + 2*ticks_interval
        ax.yaxis.set_major_locator(ticker.FixedLocator(np.arange(0,max_pos,ticks_interval)))
        r_labels = [item.get_text() for item in ax.get_yticklabels()]
        for i in range(len(r_labels)):
            if i % 2: r_labels[i]=''
        ax.set_yticklabels(r_labels, verticalalignment = "top")

        for i,lab in enumerate(fluxes):
            ax.text(th_LO[i],rr_LO[i],str(int(lab)),color='black',fontsize=11)

        # Legend
        ax.legend(loc='lower left', bbox_to_anchor=(1, 0))
        plt.show()
