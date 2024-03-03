import matplotlib
from string import ascii_uppercase
from .baseSimulation import *
import math
from scipy.optimize import curve_fit, minimize

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
                field_simul_data = np.load(file_field_simul, allow_pickle=True)
                self.asterismsRecArray = field_simul_data
                print('Number of Fields:', len(self.asterismsRecArray))
                # self.asterismsRecArrayKeys = [x[0][1] for x in self.asterismsRecArray[0].dtype.descr]
                isMono = self.asterismMode[4:]=='Mono'
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

        print('total good asterisms: ', self.cumAstSizes)

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


    def plotFieldInterval(self, fieldIndex1, fieldIndex2):
        fieldsize = 0
        fieldIndices = list(range(fieldIndex1, fieldIndex2, 1))
        for ii in fieldIndices:
            fieldsize += self.nfieldsSizes[ii]
        print("*")
        print("* Plotting Fields " + str(fieldIndex1) + ' to ' + str(fieldIndex2) + " - number of asterisms :" + str(fieldsize))
        print("*")
        self.getSourcesData(fieldIndices)
        base = self.cumAstSizes[fieldIndex1]
        covsarray = np.array(self.cov_ellipses_Asterism)[base:base+fieldsize, :,1]**2 + np.array(self.cov_ellipses_Asterism)[base:base+fieldsize, :,2]**2
        print('covsarray.shape', covsarray.shape)
        print('self.cov_ellipses_Asterism.shape', self.cov_ellipses_Asterism.shape)
        
        if covsarray.shape[0]!=0:
            dd = np.sqrt(np.average(covsarray, axis=1))
            self.twoPlots(dd, base)


    def plotField(self, fieldIndex):
        fieldsize = self.nfieldsSizes[fieldIndex]
        print("*")
        print("* Plotting Field " + str(fieldIndex) + " - number of asterisms :" + str(fieldsize))
        print("*")
        self.getSourcesData([fieldIndex])
        base = self.cumAstSizes[fieldIndex]
        covsarray = np.array(self.cov_ellipses_Asterism)[base:base+fieldsize, :,1]**2 + np.array(self.cov_ellipses_Asterism)[base:base+fieldsize, :,2]**2
        print('covsarray.shape', covsarray.shape)
        print('self.cov_ellipses_Asterism.shape', self.cov_ellipses_Asterism.shape)

        if covsarray.shape[0]!=0:
            dd = np.sqrt(np.average(np.abs(covsarray), axis=1))
            self.twoPlots(dd, base)
        
#            if not isinstance(self.asterismsRecArray[field], np.recarray):
#                continue
#            fieldsize = min( len(self.currentFieldAsterismsIndices), self.asterismsRecArray[field].size)
#                if not isinstance(self.asterismsRecArray[field][ast], np.record):
#                    continue
#                    print(self.asterismsRecArray[field][ast])
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


    def plotAsterisms(self, dd, ind_sort, istart=None, istop=None, min_id=None):
        al = np.log(dd)
        print('np.max(al)', np.max(al))
        al = (np.max(al) - al)/np.max(al)
        np.random.seed(12345)
        if not istop:
            istop = self.asterismsInputDataCartesian.shape[0]
        if not istart:
            istart = 0
        X = self.asterismsInputDataCartesian
        xcoords = X[istart:istop, 0, :]
        ycoords = X[istart:istop, 1, :]
        fluxes  = X[istart:istop, 2, :]
        max_flux = np.max(fluxes)
        print('Max flux star value: ', max_flux )
        print('min_id', min_id)
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
        
        for i in range(istop-1, istart-1, -1):
            jj = ind_sort[i-istart]
            aa = 1.0 # al[i-istart] * 0.9 + 0.1
            coords = np.transpose(X[jj+istart, :2, :])
            xx = np.sum(X[jj+istart, 0, :])/float(X.shape[2])
            yy = np.sum(X[jj+istart, 1, :])/float(X.shape[2])
            px = X[jj+istart, 0, :]
            py = X[jj+istart, 1, :]
            if X.shape[2]==1:
                circle1 = plt.Circle((xx, yy), scales[jj], color=cm1(norm(al[jj])), fill=False)
                ax.add_patch(circle1)
            else:
                ppp=ax.quiver([xx, xx, xx], [yy, yy, yy], px-xx, py-yy, color=cm1(norm(al[jj])), width=0.003, scale_units='xy', scale=1, alpha = aa )

        # draw best asterism
        coords = np.transpose(X[min_id, :2, :])
        if X.shape[2]==1:
            circle1 = plt.Circle((coords[0,0], coords[0,1]), scales[min_id-istart], color='r', fill=True, alpha=0.5)
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

        
    def computeHeuristic(self,  dd, ind_sort, istart=None, istop=None, min_id=None):

        def funcPolar(X, A, B, C, D, E, F, G, H, I, J):
            r, x, y, f = X
            #f = np.log(1/af)
            return (A * r ** 4) + (B * f ** 4) + (C * r ** 3) + (D * f ** 3) + (E * r ** 2) + (F * f ** 2) + (G * r) + (H * f) + I + (J * r * f)

        def costFunction(params, x, y):
            A, B, C, D, E, F, G, H, I, J = params
            return np.sum((y - funcPolar(x, A, B, C, D, E, F, G, H, I, J))**2)

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

        if not istop:
            istop = self.asterismsInputDataCartesian.shape[0]
        if not istart:
            istart = 0
        print('istart, istop: ', istart, istop)

        Z = np.log(dd)
        fluxes  = -np.log(self.asterismsInputDataCartesian[istart:istop, 2, :])
        rcoords = self.asterismsInputDataPolar[istart:istop, 0, :]
        X = self.asterismsInputDataCartesian
        xcoords = X[istart:istop, 0, :]
        ycoords = X[istart:istop, 1, :]
        realDataZ = Z

        func = funcMix
        realDataX = np.array( [rcoords[:,0], xcoords[:,0], ycoords[:,0], fluxes[:,0]] )

        for func in [funcPolar]:
#        for func in [funcPolar, funcCartesian, funcMix]:
            popt, pcov = curve_fit(func, realDataX, realDataZ)
            print(popt, pcov)
            for i, j in zip(popt, ascii_uppercase):
                print(f"{j} = {i:.6f}")
            zApprox = func(realDataX, *popt)
            relativeAbsoluteError = np.abs((realDataZ-zApprox)/Z)
            print( "Average Relative Absolute Error", np.mean(relativeAbsoluteError))
            plt.scatter(rcoords[:,0], relativeAbsoluteError)
            plt.show()
            plt.scatter(fluxes[:,0], relativeAbsoluteError)
            plt.show()

            initial_guess = (1,1,1,1,1,1,1,1,1,1)

            # Perform the minimization
            result = minimize(costFunction, initial_guess, args=(realDataX, realDataZ), method='Nelder-Mead' )
            print(result)
            popt = result.x
            zApprox = func(realDataX, *popt)
            relativeAbsoluteError = np.abs((realDataZ-zApprox)/Z)
            print( "Average Relative Absolute Error", np.mean(relativeAbsoluteError))
            plt.scatter(rcoords[:,0], relativeAbsoluteError)
            plt.show()
            plt.scatter(fluxes[:,0], relativeAbsoluteError)
            plt.show()


    def twoPlots(self, dd0, base):
        min_id = np.argmin(dd0)
        ind_sort = np.argsort(dd0, axis=0)
        RR = self.asterismsInputDataPolar[ind_sort, 0, 0]
        FF = self.asterismsInputDataPolar[ind_sort, 2, 0]
        # print(dd0)
        # print(ind_sort)
        # print(dd0[ind_sort])       
        # print(RR)
        # print(FF)
        plt.plot(dd0[ind_sort])
        # plt.plot(RR)
        # plt.plot(FF)
        # plt.gca().set_yscale("log")
        # plt.gca().set_yscale("linear")
        plt.show()
        #print('base, min_id ', base, min_id)
        self.plotAsterisms(dd0, ind_sort, base, base+dd0.shape[0], base+min_id)
        self.computeHeuristic(dd0, ind_sort, base, base+dd0.shape[0], base+min_id)

