import matplotlib
from string import ascii_uppercase
from .baseSimulation import *
import math
from scipy.optimize import curve_fit, minimize


def funcPolar(X, A0, A, B, C, D, E0, E, F, G, H, I, J, J0, J1):
    r, x, y, af = X
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


class asterismSimulation(baseSimulation):


    def __init__(self, simulName, path, parametersFile, outputDir,
                 outputFile, doPlot=False, addSrAndFwhm=False, verbose=False,
                 getHoErrorBreakDown=False):

        super().__init__(path, parametersFile, outputDir, outputFile, doConvolve=True,
                          doPlot=False, addSrAndFwhm=addSrAndFwhm,
                          verbose=verbose, getHoErrorBreakDown=False,
                          savePSDs=False)

        self.simulName = simulName
        # store the parameters in data members
        self.asterismsInputDataCartesian = None
        self.asterismsInputDataPolar = None
        self.hasAsterismSection = False
        if 'ASTERISM_SELECTION' in self.my_data_map.keys():
            self.hasAsterismSection = True
            self.globalOffset = 0
            self.asterismMode = self.my_data_map['ASTERISM_SELECTION']['mode']
            listZ = self.my_data_map['ASTERISM_SELECTION']['Zenith']
            listA = self.my_data_map['ASTERISM_SELECTION']['Azimuth']
            listP = self.my_data_map['ASTERISM_SELECTION']['NumberPhotons']
            listF  = self.my_data_map['ASTERISM_SELECTION']['Frequencies']
            if 'heuristicModel' in self.my_data_map['ASTERISM_SELECTION'].keys():
                heuristicModel = self.my_data_map['ASTERISM_SELECTION']['heuristicModel']
            else:
                heuristicModel = None
            if not isinstance(listF, list):
                listF  = [listF] * len(listZ)
            self.cumAstSizes = [0]
            self.nfields = 1
            if self.asterismMode=='Sets':
                pointings = polarToCartesian(np.array( [listZ, listA]))
                xxPointigs  = pointings[0,:]
                yyPointigs  = pointings[1,:]
                self.nfieldsSizes = [pointings.shape[1]]
                self.asterismsInputDataCartesian = np.asarray( [xxPointigs, yyPointigs, listP, listF ] )
                self.asterismsInputDataCartesian = np.swapaxes(self.asterismsInputDataCartesian, 0,1)                
                self.asterismsInputDataPolar = np.asarray( [listZ, listA, listP, listF ] )
                self.asterismsInputDataPolar = np.swapaxes(self.asterismsInputDataPolar, 0,1)
            elif self.asterismMode[:7]=='Singles':
                if self.asterismMode[7]=='3' or self.asterismMode[7]=='1':
                    nStars = len(self.my_data_map['ASTERISM_SELECTION']['Zenith'])
                    pointings = polarToCartesian(np.array( [listZ, listA]))
                    xxPointigs  = pointings[0,:]
                    yyPointigs  = pointings[1,:]
                else:
                    self.asterismMode = 'INVALID'
                    self.hasAsterismSection = False
                    print('ERROR: Only Singles1 (One Start Asterisms) and Singles3 (3 Stars Asterisms) are implemented.')
                    return

                all_combos = list(itertools.combinations(list(range(nStars)), int(self.asterismMode[7])))
                self.nfieldsSizes = [len(all_combos)]
                self.cumAstSizes.append(self.nfieldsSizes[0])
                self.asterismsInputDataCartesian = np.array( [np.take(xxPointigs, all_combos),
                                                     np.take(yyPointigs, all_combos),
                                                     np.take(np.array(listP, dtype=np.float64), all_combos),
                                                     np.take(np.array(listF, dtype=np.float64), all_combos)] )
                self.asterismsInputDataCartesian = np.swapaxes(self.asterismsInputDataCartesian, 0,1) 
                self.asterismsInputDataPolar = np.array( [np.take(np.array(listZ, dtype=np.float64), all_combos),
                                                     np.take(np.array(listA, dtype=np.float64), all_combos),
                                                     np.take(np.array(listP, dtype=np.float64), all_combos),
                                                     np.take(np.array(listF, dtype=np.float64), all_combos)] )
                self.asterismsInputDataPolar = np.swapaxes(self.asterismsInputDataPolar, 0,1)
                self.currentFieldAsterismsIndices = all_combos
            elif self.asterismMode=='Generate':
                self.asterismsInputDataCartesian, self.asterismsInputDataPolar = self.generateTriangles(listZ[0], listZ[0], 10)
                
            elif self.asterismMode[:4]=='File':
                self.nfields = self.my_data_map['ASTERISM_SELECTION']['fieldsNumber']
                # number of asterisms for each field
                self.nfieldsSizes = []
                file_field_simul = self.my_data_map['ASTERISM_SELECTION']['filename']
                self.globalOffset = self.my_data_map['ASTERISM_SELECTION']['offset']
                # self.asterismsRecArrayKeys = [x[0][1] for x in self.asterismsRecArray[0].dtype.descr]
                isMono = self.asterismMode[-4:]=='Mono'
                if self.asterismMode[4:10]=='Random':
                    if file_field_simul=='':
                        self.asterismsInputDataCartesian, self.asterismsInputDataPolar = self.generateRandom(self.nfields, isMono, 'savedRandom')
                    else:
                        self.asterismsInputDataCartesian = np.load(os.path.join(self.outputDir, file_field_simul +'C.npy'))
                        self.asterismsInputDataPolar = np.load(os.path.join(self.outputDir, file_field_simul +'P.npy'))
                        self.nfieldsSizes = np.load(os.path.join(self.outputDir, file_field_simul +'F.npy')).tolist()
                        self.cumAstSizes = np.load(os.path.join(self.outputDir, file_field_simul +'S.npy')).tolist()
                else:
                    field_simul_data = np.load(file_field_simul, allow_pickle=True)
                    self.asterismsRecArray = field_simul_data
                    print('Number of Fields:', len(self.asterismsRecArray))
                    self.asterismsInputDataCartesian, self.asterismsInputDataPolar = self.generateFromRecArray(self.nfields, isMono)
            else:
                self.asterismMode = 'INVALID'                
        else:
            self.hasAsterismSection = False


    def configLO(self, astIndex=None):
        if astIndex==0:
            self.cartSciencePointingCoords = np.dstack( (self.xxSciencePointigs, self.yySciencePointigs) ).reshape(-1, 2)
            # Here we assume the same wavelenght for all the phon counts of the stars in the asterism
            LO_wvl_temp = self.my_data_map['sources_LO']['Wavelength']
            if isinstance(LO_wvl_temp, list):
                self.LO_wvl = LO_wvl_temp[0]  # lambda
            else:
                self.LO_wvl = LO_wvl_temp     # lambda
            self.my_data_map['sources_LO']['Zenith'] = self.currentFieldSourcesData['Zenith']
            self.my_data_map['sources_LO']['Azimuth'] = self.currentFieldSourcesData['Azimuth']
            self.my_data_map['sensor_LO']['NumberPhotons'] = self.currentFieldSourcesData['NumberPhotons']                
            self.my_data_map['RTC']['SensorFrameRate_LO'] = self.currentFieldSourcesData['Frequencies']
            self.LO_zen_field    = self.my_data_map['sources_LO']['Zenith']
            self.LO_az_field     = self.my_data_map['sources_LO']['Azimuth']
            self.LO_fluxes_field = self.my_data_map['sensor_LO']['NumberPhotons']
            self.LO_freqs_field = self.my_data_map['RTC']['SensorFrameRate_LO']
            self.NGS_fluxes_field = []
            polarNGSCoordsList = []
            for aFr, aFlux, aZen, aAz in zip(self.LO_freqs_field, self.LO_fluxes_field, self.LO_zen_field, self.LO_az_field):
                polarNGSCoordsList.append([aZen, aAz])
                self.NGS_fluxes_field.append(aFlux*aFr)
            polarNGSCoords     = np.asarray(polarNGSCoordsList)
            self.nNaturalGS_field  = len(self.LO_zen_field)
            cartNGSCoordsList = []
            for i in range(self.nNaturalGS_field):
                cartNGSCoordsList.append(polarToCartesian(polarNGSCoords[i,:]))
            self.cartNGSCoords_field = np.asarray(cartNGSCoordsList)

        self.currentAsterismIndices = self.currentFieldAsterismsIndices[astIndex]
        super().setAsterismData()
       
            
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


    def freqFromMagnitudeERIS(self, magnitude):
        if magnitude>17:
            return 100.0
        elif magnitude>16:
            return 250.0
        else:
            return 500.0

    def freqsFromMagnitudes(self, magnitudes):
        result = []
        for m in magnitudes:
            result.append(self.freqFromMagnitudeERIS(m))
        return result

    def fluxFromMagnitudeERIS(self, m):
        m0 = 0.0
        f0 = 1.51e10
        return f0 * np.power(10.0, (-(m-m0)/2.5) )


    def generateRandom(self, max_field, isMono, filename):
        self.setsList = []
        self.polarSetsList = []
        self.nfieldsSizes = []
        self.cumAstSizes = [0]
        self.skippedFieldIndexes = []
        total_skipped_fields = 0
        total_skipped_asterisms = 0
        ObscurationRatio   = self.my_data_map['telescope']['ObscurationRatio']
        TelescopeDiameter  = self.my_data_map['telescope']['TelescopeDiameter']
        NumberLenslets     = self.my_data_map['sensor_LO']['NumberLenslets']
        transmissionFactor = self.my_data_map['ASTERISM_SELECTION']['transmissionFactor']
        bands              = self.my_data_map['ASTERISM_SELECTION']['bands']
        freqs                = self.my_data_map['ASTERISM_SELECTION']['Frequencies']
        N_sa_tot_LO        = NumberLenslets[0]**2
        if NumberLenslets[0] > 2:
            N_sa_tot_LO   = int ( np.floor( N_sa_tot_LO * np.pi/4.0 * (1.0 - ObscurationRatio**2) ) )

        for i in range(self.globalOffset, self.globalOffset+max_field, 1):
            asterismSet = {}
            number_of_asterisms0 = np.random.randint(1,20) # len(self.asterismsRecArray[i])
            number_of_asterisms = number_of_asterisms0
            if isMono:
                number_of_asterisms = 0
                self.currentFieldSourcesData = {}
                self.currentFieldSourcesData['Zenith'] = []
                self.currentFieldSourcesData['Azimuth'] = []
                self.currentFieldSourcesData['NumberPhotons'] = []
                self.currentFieldSourcesData['Frequencies'] = []
            for j in range(number_of_asterisms0):
                # asterism=None
                xcoords = np.random.uniform(-30,30) # self.asterismsRecArray[i][j]['COORD'][0]
                ycoords = np.random.uniform(-30,30) # self.asterismsRecArray[i][j]['COORD'][1]
                flux = 0.0
                magnitude = 0.0
                for b in bands:
                    mm = np.random.uniform(11,20)
                    magnitude += mm
                    flux += self.fluxFromMagnitudeERIS(mm)

                magnitude /= float(len(bands))
                freqs = self.freqsFromMagnitudes([magnitude])[0]
                fluxScaling = ((TelescopeDiameter/2.0)**2 * np.pi * (1-ObscurationRatio**2) * transmissionFactor / N_sa_tot_LO)
                flux = flux * fluxScaling / freqs
                if np.min(flux)<=0.0:
#                            print("Field:" + str(i) + "- Asterism:" + str(j) + " SKIPPED because of Flux 0 star")
                    total_skipped_asterisms += 1
                    if not isMono:
                        number_of_asterisms -= 1
                    # self.skippedAsterismsIndexes.append(j+self.cumAstSizes[-1])
                    continue
                if isMono:
                    source = np.array([xcoords, ycoords, flux, freqs])
                    self.currentFieldSourcesData['Zenith'].append(source[0])
                    self.currentFieldSourcesData['Azimuth'].append(source[1])
                    self.currentFieldSourcesData['NumberPhotons'].append(source[2])
                    self.currentFieldSourcesData['Frequencies'].append(source[3])
                    asterism = np.vstack( [[xcoords], [ycoords], [flux], [freqs]] )
                    pasterism = np.vstack( [cartesianToPolar(asterism[0:2,:]), [flux], [freqs] ])
                    self.setsList.append(asterism)
                    self.polarSetsList.append(pasterism)
                    number_of_asterisms += 1
                else:
                    asterism = np.vstack( [xcoords, ycoords, flux, freqs] )
                    pasterism = np.vstack( (cartesianToPolar(asterism[0:2,:]), flux, freqs) )
#                        if not asterism is None:
                    self.setsList.append(asterism)
                    self.polarSetsList.append(pasterism)

            number_of_asterisms = max(0, number_of_asterisms)
            self.nfieldsSizes.append(number_of_asterisms)
            self.cumAstSizes.append(self.cumAstSizes[-1]+number_of_asterisms)
        print('total_skipped_fields: ', total_skipped_fields)
        print('total_skipped_asterisms: ', total_skipped_asterisms)
        print('total good asterisms: ', self.cumAstSizes[-1])
        print('total good asterisms: ', self.cumAstSizes)
        np.save(os.path.join(self.outputDir, filename +'C.npy'), np.array(self.setsList))
        np.save(os.path.join(self.outputDir, filename +'P.npy'), np.array(self.polarSetsList))
        np.save(os.path.join(self.outputDir, filename +'F.npy'), np.array(self.nfieldsSizes))
        np.save(os.path.join(self.outputDir, filename +'S.npy'), np.array(self.cumAstSizes))
        return np.array(self.setsList), np.array(self.polarSetsList)


    def generateFromRecArray(self, max_field, isMono):
        self.setsList = []
        self.polarSetsList = []
        self.nfieldsSizes = []
        self.cumAstSizes = [0]
        self.skippedFieldIndexes = []
#        self.skippedAsterismsIndexes = []
        total_skipped_fields = 0
        total_skipped_asterisms = 0
        ObscurationRatio   = self.my_data_map['telescope']['ObscurationRatio']
        TelescopeDiameter  = self.my_data_map['telescope']['TelescopeDiameter']
        NumberLenslets     = self.my_data_map['sensor_LO']['NumberLenslets']
        transmissionFactor = self.my_data_map['ASTERISM_SELECTION']['transmissionFactor']
        bands              = self.my_data_map['ASTERISM_SELECTION']['bands']
        freqs                = self.my_data_map['ASTERISM_SELECTION']['Frequencies']
        N_sa_tot_LO        = NumberLenslets[0]**2
        if NumberLenslets[0] > 2:
            N_sa_tot_LO   = int ( np.floor( N_sa_tot_LO * np.pi/4.0 * (1.0 - ObscurationRatio**2) ) )

        for i in range(self.globalOffset, self.globalOffset+max_field, 1):
            skipped_field = False
            if type(self.asterismsRecArray[i]) is np.int16 or type(self.asterismsRecArray[i]) is np.int64:
                print("Field:" + str(i) + " SKIPPED")
                total_skipped_fields += 1
                self.nfields -=1
                skipped_field = True
                self.skippedFieldIndexes.append(i)
            else:
                if skipped_field:
                    number_of_asterisms0 = 0
                else:
                    asterismSet = {}
                    number_of_asterisms0 = len(self.asterismsRecArray[i])
                    number_of_asterisms = number_of_asterisms0
                    if isMono:
                        number_of_asterisms = 0
                        self.currentFieldSourcesData = {}
                        self.currentFieldSourcesData['Zenith'] = []
                        self.currentFieldSourcesData['Azimuth'] = []
                        self.currentFieldSourcesData['NumberPhotons'] = []
                        self.currentFieldSourcesData['Frequencies'] = []
                    for j in range(number_of_asterisms0):
                        # asterism=None
                        xcoords = self.asterismsRecArray[i][j]['COORD'][0]
                        ycoords = self.asterismsRecArray[i][j]['COORD'][1]
                        flux = 0.0
                        magnitude = 0.0
                        for b in bands:
                            flux += self.asterismsRecArray[i][j]['FLUX' + b]
                            magnitude += self.asterismsRecArray[i][j][b+'MAG']
                        magnitude /= float(len(bands))
                        freqs = self.freqsFromMagnitudes(magnitude)
                        fluxScaling = ((TelescopeDiameter/2.0)**2 * np.pi * (1-ObscurationRatio**2) * transmissionFactor / N_sa_tot_LO)
                        flux = flux * fluxScaling / freqs
                        if np.min(flux)<=0.0:
#                            print("Field:" + str(i) + "- Asterism:" + str(j) + " SKIPPED because of Flux 0 star")
                            total_skipped_asterisms += 1
                            if not isMono:
                                number_of_asterisms -= 1
                            # self.skippedAsterismsIndexes.append(j+self.cumAstSizes[-1])
                            continue
                        if isMono:
                            for si in range(3):
                                source = np.array([xcoords[si], ycoords[si], flux[si], freqs[si]])
                                s_index = self.sourceIsPresent(source)
                                if s_index==-1 and np.abs(source[0])<30.0 and np.abs(source[1])<30.0:
                                    self.currentFieldSourcesData['Zenith'].append(source[0])
                                    self.currentFieldSourcesData['Azimuth'].append(source[1])
                                    self.currentFieldSourcesData['NumberPhotons'].append(source[2])
                                    self.currentFieldSourcesData['Frequencies'].append(source[3])
#                                    if magnitude[si]>19.0 or flux[si]<=0.0:
#                                        total_skipped_asterisms += 1
#                                        #number_of_asterisms -= 1
#                                    else:

                                    asterism = np.vstack( [[xcoords[si]], [ycoords[si]], [flux[si]], [freqs[si]]] )
                                    pasterism = np.vstack( [cartesianToPolar(asterism[0:2,:]), [flux[si]], [freqs[si]] ])
                                    self.setsList.append(asterism)
                                    self.polarSetsList.append(pasterism)
                                    number_of_asterisms += 1
                        else:
                            asterism = np.vstack( [xcoords, ycoords, flux, freqs] )
                            pasterism = np.vstack( (cartesianToPolar(asterism[0:2,:]), flux, freqs) )
#                        if not asterism is None:
                            self.setsList.append(asterism)
                            self.polarSetsList.append(pasterism)

                    number_of_asterisms = max(0, number_of_asterisms)
                    self.nfieldsSizes.append(number_of_asterisms)
                    self.cumAstSizes.append(self.cumAstSizes[-1]+number_of_asterisms)

        print('total_skipped_fields: ', total_skipped_fields)
        print('total_skipped_asterisms: ', total_skipped_asterisms)
        print('total good asterisms: ', self.cumAstSizes[-1])

        return np.array(self.setsList), np.array(self.polarSetsList)


    def getAsterismSource(self, field, ast, s):
        astIndexGlobal = self.cumAstSizes[field]+ast
        z = self.asterismsInputDataPolar[astIndexGlobal, 0, s]
        a = self.asterismsInputDataPolar[astIndexGlobal, 1, s]
        n = self.asterismsInputDataPolar[astIndexGlobal, 2, s]
        f = self.asterismsInputDataPolar[astIndexGlobal, 3, s]
        return np.array([z, a, n, f])


    def sourceIsPresent(self, source):
        # a source is a triple (array) of zenith, azimuth, photons
        for z, a, n, f, ii in zip(self.currentFieldSourcesData['Zenith'],
                               self.currentFieldSourcesData['Azimuth'], 
                               self.currentFieldSourcesData['NumberPhotons'],
                               self.currentFieldSourcesData['Frequencies'],
                               range(len(self.currentFieldSourcesData['Zenith'])) ):
            if np.allclose(np.array([z,a,n,f]), source):
                return ii
        return -1


    def getSourcesData(self, fields):
        self.currentFieldSourcesData = {}
        self.currentFieldSourcesData['Zenith'] = [] 
        self.currentFieldSourcesData['Azimuth'] = []
        self.currentFieldSourcesData['NumberPhotons'] = []
        self.currentFieldSourcesData['Frequencies'] = []
        self.currentFieldAsterismsIndices = []
        for field in fields:
            fieldsize = self.nfieldsSizes[field]
#            print('fieldsize', fieldsize)
            for ast in range(fieldsize):
                astIndices = []
                starsInAsterism = self.asterismsInputDataPolar.shape[2]
    #            print('starsInAsterism', starsInAsterism)
                for s in range(starsInAsterism):
                    source = self.getAsterismSource(field, ast, s)
                    s_index = self.sourceIsPresent(source)
                    if s_index==-1:
                        self.currentFieldSourcesData['Zenith'].append(source[0])
                        self.currentFieldSourcesData['Azimuth'].append(source[1])
                        self.currentFieldSourcesData['NumberPhotons'].append(source[2])
                        self.currentFieldSourcesData['Frequencies'].append(source[3])
                        s_index = len(self.currentFieldSourcesData['Zenith'])-1
                    astIndices.append(s_index)
                self.currentFieldAsterismsIndices.append(astIndices)


    def selectData(self, fieldIndex1, fieldIndex2 = None):
        if fieldIndex2 is None:
            self.currentFieldsize = self.nfieldsSizes[fieldIndex1]
            self.getSourcesData([fieldIndex1])
            self.currentBase = self.cumAstSizes[fieldIndex1]
            self.covsarray = np.array(self.cov_ellipses_Asterism)[self.currentBase:self.currentBase+self.currentFieldsize, :,1]**2 + np.array(self.cov_ellipses_Asterism)[self.currentBase:self.currentBase+self.currentFieldsize, :,2]**2
        else:
            self.currentFieldsize = 0
            fieldIndices = list(range(fieldIndex1, fieldIndex2, 1))
            for ii in fieldIndices:
                self.currentFieldsize += self.nfieldsSizes[ii]
            self.getSourcesData(fieldIndices)
            self.currentBase = self.cumAstSizes[fieldIndex1]
            self.covsarray = np.array(self.cov_ellipses_Asterism)[self.currentBase:self.currentBase+self.currentFieldsize, :,1]**2 + np.array(self.cov_ellipses_Asterism)[self.currentBase:self.currentBase+self.currentFieldsize, :,2]**2
        if self.covsarray.shape[0]!=0:
            self.jitterMeasure = np.sqrt(np.average(self.covsarray, axis=1))
            self.minJitter_id = np.argmin(self.jitterMeasure)+self.currentBase
            self.sortedJitterIndices = np.argsort(self.jitterMeasure, axis=0)
            self.lastJitterIndex = self.jitterMeasure.shape[0] + self.currentBase


    def fitHeuristicModel(self, fieldIndex1, fieldIndex2, modelName):
        self.selectData(fieldIndex1, fieldIndex2)
        self.computeHeuristicModel(modelName)


    def testHeuristicModel(self, fieldIndex1, fieldIndex2, modelName):
        self.selectData(fieldIndex1, fieldIndex2)
        popt = np.load(os.path.join(self.outputDir, modelName +'.npy'))
        Z = np.log(self.jitterMeasure+1)
        # Z = self.jitterMeasure
        rcoords = self.asterismsInputDataPolar[self.currentBase:self.lastJitterIndex, 0, :]
        freqs   = self.asterismsInputDataCartesian[self.currentBase:self.lastJitterIndex, 3, :]
        fluxes  = self.asterismsInputDataCartesian[self.currentBase:self.lastJitterIndex, 2, :]
        fluxes  = np.log(fluxes*freqs+1.0)
        X = self.asterismsInputDataCartesian
        xcoords = X[self.currentBase:self.lastJitterIndex, 0, :]
        ycoords = X[self.currentBase:self.lastJitterIndex, 1, :]
        realDataZ = Z
        func = funcPolar
        realDataX = np.array( [rcoords[:,0], xcoords[:,0], ycoords[:,0], fluxes[:,0]] )
        # popt, pcov = curve_fit(func, realDataX, realDataZ)
#        for i, j in zip(popt, ascii_uppercase):
#            print(f"{j} = {i:.6f}")
        zApprox = func(realDataX, *popt)
        relativeAbsoluteError = (realDataZ-zApprox)/Z
        sortedJitterIndicesModel = np.argsort(zApprox, axis=0)
        print( "Average Relative Absolute Error", np.mean(np.abs(relativeAbsoluteError)))
        plt.scatter(rcoords[:,0], relativeAbsoluteError)
        plt.show()
        plt.scatter(fluxes[:,0], relativeAbsoluteError)
        plt.show()
        plt.scatter(Z, relativeAbsoluteError)
        plt.show()
        totalS = 0
        for field in range(fieldIndex1, fieldIndex2, 1):
            self.currentField = field
            self.getSourcesData([field])
            fieldsize = len(self.currentFieldAsterismsIndices)
            base = self.cumAstSizes[field]

            self.selectData(field)
            Z = np.log(np.abs(self.jitterMeasure)+1.0)
            rcoords = self.asterismsInputDataPolar[self.currentBase:self.lastJitterIndex, 0, :]
            freqs   = self.asterismsInputDataCartesian[self.currentBase:self.lastJitterIndex, 3, :]
            fluxes  = self.asterismsInputDataCartesian[self.currentBase:self.lastJitterIndex, 2, :]
            fluxes  = np.log(fluxes*freqs+1.0)
            X = self.asterismsInputDataCartesian
            xcoords = X[self.currentBase:self.lastJitterIndex, 0, :]
            ycoords = X[self.currentBase:self.lastJitterIndex, 1, :]
            realDataZ = Z
            func = funcPolar
            realDataX = np.array( [rcoords[:,0], xcoords[:,0], ycoords[:,0], fluxes[:,0]] )
            zApprox = func(realDataX, *popt)
            if zApprox.shape[0]>1:
                relativeAbsoluteError = (realDataZ-zApprox)/Z
                sortedJitterIndicesModel = np.argsort(zApprox, axis=0)
                zApproxSorted = zApprox[sortedJitterIndicesModel]
                deltas = np.abs( 2.0 * (zApproxSorted[:-1] - zApproxSorted[1:]) / ( zApproxSorted[1:] + zApproxSorted[:-1]))
                cumDeltas = np.cumsum(deltas)
                if deltas[0]<0.03:
                    ss = np.where(cumDeltas<0.04)[0].shape[0]+1
                    print(field, 'Low Delta', ss)
                    totalS += ss
                #if not np.equal(self.sortedJitterIndices, sortedJitterIndicesModel).all():
                if self.sortedJitterIndices[0] != sortedJitterIndicesModel[0]:
                    print(field, deltas)
                    print(self.sortedJitterIndices)
                    print(sortedJitterIndicesModel)
                    print(zApproxSorted)
                    if not deltas[0]<0.04:
                        print('ERROR: Correct ordering estimation missed!!!!!!!!')
        print('totalS', totalS)


    def  computeHeuristicModel(self, modelName):
        Z = np.log(self.jitterMeasure+1)
        #Z = self.jitterMeasure
        rcoords = self.asterismsInputDataPolar[self.currentBase:self.lastJitterIndex, 0, :]
        freqs  = self.asterismsInputDataCartesian[self.currentBase:self.lastJitterIndex, 3, :]
        fluxes = self.asterismsInputDataCartesian[self.currentBase:self.lastJitterIndex, 2, :]
        fluxes  = np.log(fluxes*freqs+1.0)

        X = self.asterismsInputDataCartesian
        xcoords = X[self.currentBase:self.lastJitterIndex, 0, :]
        ycoords = X[self.currentBase:self.lastJitterIndex, 1, :]
        realDataZ = Z
        func = funcPolar
        realDataX = np.array( [rcoords[:,0], xcoords[:,0], ycoords[:,0], fluxes[:,0]] )
        for func in [funcPolar]:
            popt, pcov = curve_fit(func, realDataX, realDataZ)
            print(popt, pcov)
            for i, j in zip(popt, ascii_uppercase):
                print(f"{j} = {i:.6f}")
            zApprox = func(realDataX, *popt)
            relativeAbsoluteError = (realDataZ-zApprox)/Z
            print( "Average Relative Absolute Error", np.mean(np.abs(relativeAbsoluteError)))
            plt.scatter(rcoords[:,0], relativeAbsoluteError)
            plt.show()
            plt.scatter(fluxes[:,0], relativeAbsoluteError)
            plt.show()
            plt.scatter(Z, relativeAbsoluteError)
            plt.show()
            np.save(os.path.join(self.outputDir, modelName +'.npy'), np.array(popt))


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
        print('self.cov_ellipses_Asterism.shape', self.cov_ellipses_Asterism.shape)
        if self.covsarray.shape[0]!=0:
            self.twoPlots()


    def computeAsterisms(self, eeRadiusInMas):
        self.sr_Asterism = []
        self.fwhm_Asterism = []
        self.ee_Asterism = []
        self.cov_ellipses_Asterism = []
        for field in range(self.nfields):
            self.currentField = field
            self.getSourcesData([field])
            fieldsize = len(self.currentFieldAsterismsIndices)
            base = self.cumAstSizes[field]
            for ast in range(fieldsize):
                self.currentAsterism = ast
#                try:
                self.doOverallSimulation(ast)
                self.computeMetrics(eeRadiusInMas)
#               except:
#                    print("Unexpected Error Computing (field, asterism): ", self.currentField, self.currentAsterism)
#                    print(self.currentFieldSourcesData)
                self.sr_Asterism.append(self.sr)
                self.fwhm_Asterism.append(self.fwhm) 
                self.ee_Asterism.append(self.ee) 
                self.cov_ellipses_Asterism.append(self.cov_ellipses)
            if (field+1) % 100 == 0 or field==self.nfields-1:
                print("Field " + str(field) + " DONE")
                np.save(os.path.join(self.outputDir, self.simulName+'fw.npy'), np.array(self.fwhm_Asterism))
                np.save(os.path.join(self.outputDir, self.simulName+'ee.npy'), np.array(self.ee_Asterism))
                np.save(os.path.join(self.outputDir, self.simulName+'covs.npy'), np.array(self.cov_ellipses_Asterism))
                np.save(os.path.join(self.outputDir, self.simulName+'sr.npy'), np.array(self.sr_Asterism))


    def reloadResults(self):
        self.fwhm_Asterism = np.load(os.path.join(self.outputDir, self.simulName+'fw.npy'))
        self.ee_Asterism = np.load(os.path.join(self.outputDir, self.simulName+'ee.npy'))
        self.cov_ellipses_Asterism = np.load(os.path.join(self.outputDir, self.simulName+'covs.npy'))
        self.sr_Asterism = np.load(os.path.join(self.outputDir, self.simulName+'sr.npy'))


    def plotAsterisms(self):
        al = np.log(self.jitterMeasure)
        print('np.max(al)', np.max(al))
        al = (np.max(al) - al)/np.max(al)
        np.random.seed(12345)
        if not self.currentBase:
            self.currentBase = 0
        X = self.asterismsInputDataCartesian
        xcoords = X[self.currentBase:self.lastJitterIndex, 0, :]
        ycoords = X[self.currentBase:self.lastJitterIndex, 1, :]
        fluxes  = X[self.currentBase:self.lastJitterIndex, 2, :]
        max_flux = np.max(fluxes)
        print('Max flux star value: ', max_flux )
        print('self.minJitter_id', self.minJitter_id)
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
        ax.set_xlim([-60,60]) 
        ax.set_ylim([-60,60]) 
        cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
        cbar = fig.colorbar(sm, cax=cax) # Similar to fig.colorbar(im, cax = cax)
        cbar.set_label('Asterim Rating')
        plt.show()


    def twoPlots(self):
        # RR = self.asterismsInputDataPolar[self.sortedJitterIndices, 0, 0]
        # FF = self.asterismsInputDataPolar[self.sortedJitterIndices, 2, 0]
        # print(self.jitterMeasure)
        # print(self.sortedJitterIndices)
        # print(self.jitterMeasure[self.sortedJitterIndices])
        # print(RR)
        # print(FF)
        plt.plot(self.jitterMeasure[self.sortedJitterIndices])
        # plt.plot(RR)
        # plt.plot(FF)
        # plt.gca().set_yscale("log")
        # plt.gca().set_yscale("linear")
        plt.show()
        #print('self.currentBase, self.minJitter_id ', self.currentBase, self.minJitter_id)
        self.plotAsterisms()

'''
            initial_guess = (0,0,0,0, 0,0,0,0, 0,0)
            # Perform the minimization
            result = minimize(costFunction, initial_guess, args=(realDataX, realDataZ), method='L-BFGS-B' )
            print(result)
            popt = result.x
            zApprox = func(realDataX, *popt)
            relativeAbsoluteError = (realDataZ-zApprox)/Z
            print( "Average Relative Absolute Error", np.mean(np.abs(relativeAbsoluteError)))
            plt.scatter(rcoords[:,0], relativeAbsoluteError)
            plt.show()
            plt.scatter(fluxes[:,0], relativeAbsoluteError)
            plt.show()
            plt.scatter(Z, relativeAbsoluteError)
            plt.show()
'''
