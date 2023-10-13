from .baseSimulation import *

class asterismSimulation(baseSimulation):
    
        
    def __init__(self, path, parametersFile, outputDir, outputFile,
                          doPlot=False, addSrAndFwhm=False,
                          verbose=False, getHoErrorBreakDown=False):

        super().__init__(path, parametersFile, outputDir, outputFile, doConvolve=True,
                          doPlot=False, addSrAndFwhm=addSrAndFwhm,
                          verbose=verbose, getHoErrorBreakDown=False,
                          savePSDs=False)

        # store the parameters in data members
        self.asterismsInputDataCartesian = None       
        self.asterismsInputDataPolar = None       
        self.hasAsterismSection = False
        if 'ASTERISM_SELECTION' in self.my_data_map.keys():
            self.hasAsterismSection = True
            self.asterismMode = self.my_data_map['ASTERISM_SELECTION']['mode']
            listZ = self.my_data_map['ASTERISM_SELECTION']['Zenith']
            listA = self.my_data_map['ASTERISM_SELECTION']['Azimuth']
            listP = self.my_data_map['ASTERISM_SELECTION']['NumberPhotons']
            self.cumAstSizes = [0]
            self.nfields = 1
            if self.asterismMode=='Triplets':
                pointings = polarToCartesian(np.array( [listZ, listA]))
                xxPointigs  = pointings[0,:]
                yyPointigs  = pointings[1,:]                
                self.asterismsInputDataCartesian = np.asarray( [xxPointigs, yyPointigs, listP ] )
                self.asterismsInputDataCartesian = np.swapaxes(self.asterismsInputDataCartesian, 0,1)                
                self.asterismsInputDataPolar = np.asarray( [listZ, listA, listP ] )
                self.asterismsInputDataPolar = np.swapaxes(self.asterismsInputDataPolar, 0,1)
            elif self.asterismMode=='Singles':
                nStars = len(self.my_data_map['ASTERISM_SELECTION']['Zenith'])
                pointings = polarToCartesian(np.array( [listZ, listA]))
                xxPointigs  = pointings[0,:]
                yyPointigs  = pointings[1,:]
                all_combos = list(itertools.combinations(list(range(nStars)), 3))
                self.asterismsInputDataCartesian = np.array( [np.take(xxPointigs, all_combos),
                                                     np.take(yyPointigs, all_combos),
                                                     np.take(np.array(listP, dtype=np.float64), all_combos)] )
                self.asterismsInputDataCartesian = np.swapaxes(self.asterismsInputDataCartesian, 0,1) 
                self.asterismsInputDataPolar = np.array( [np.take(np.array(listZ, dtype=np.float64), all_combos),
                                                     np.take(np.array(listA, dtype=np.float64), all_combos),
                                                     np.take(np.array(listP, dtype=np.float64), all_combos)] )
                self.asterismsInputDataPolar = np.swapaxes(self.asterismsInputDataPolar, 0,1)                
            elif self.asterismMode=='Generate':
                self.asterismsInputDataCartesian, self.asterismsInputDataPolar = self.generateTriangles(listZ[0], listZ[0], 10)
            elif self.asterismMode=='File':
                self.nfields = self.my_data_map['ASTERISM_SELECTION']['fieldsNumber']
                # number of asterisms for each field
                self.nfieldsSizes = []
                file_field_simul = self.my_data_map['ASTERISM_SELECTION']['filename']
                field_simul_data = np.load(file_field_simul, allow_pickle=True)
                self.asterismsRecArray = field_simul_data
                print('Number of Fields:', len(self.asterismsRecArray))
                self.asterismsRecArrayKeys = [x[0][1] for x in self.asterismsRecArray[0].dtype.descr]
                self.asterismsInputDataCartesian, self.asterismsInputDataPolar = self.generateFromRecArray(self.nfields)
            else:
                self.asterismMode = 'INVALID'                
        else:
            self.hasAsterismSection = False


    def configLO(self, astIndex=None):
        if astIndex==0:
            self.cartSciencePointingCoords = np.dstack( (self.xxSciencePointigs, self.yySciencePointigs) ).reshape(-1, 2)
            self.fr         = self.my_data_map['RTC']['SensorFrameRate_LO']
            # Here we assume the same wavelenght for all the phon counts of the stars in the asterism
            LO_wvl_temp = self.my_data_map['sources_LO']['Wavelength']
            if isinstance(LO_wvl_temp, list):
                self.LO_wvl = LO_wvl_temp[0]  # lambda
            else:
                self.LO_wvl = LO_wvl_temp     # lambda
            self.my_data_map['sources_LO']['Zenith'] = self.currentFieldSourcesData['Zenith']
            self.my_data_map['sources_LO']['Azimuth'] = self.currentFieldSourcesData['Azimuth']
            self.my_data_map['sensor_LO']['NumberPhotons'] = self.currentFieldSourcesData['NumberPhotons']                
            self.LO_zen_field    = self.my_data_map['sources_LO']['Zenith']
            self.LO_az_field     = self.my_data_map['sources_LO']['Azimuth']
            self.LO_fluxes_field = self.my_data_map['sensor_LO']['NumberPhotons']
            self.NGS_fluxes_field = []
            polarNGSCoordsList = []
            for aFlux, aZen, aAz in zip(self.LO_fluxes_field, self.LO_zen_field, self.LO_az_field):
                polarNGSCoordsList.append([aZen, aAz])
                self.NGS_fluxes_field.append(aFlux*self.fr)
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
        for rt in np.linspace(0, maxRadius, samples):
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
                #triangle = self.translateTriangle(triangleS, rt, standardPosRadius)
                self.trianglesList.append(triangleS)                
                pTriangle = np.vstack( (cartesianToPolar(triangleS[0:2,:]), np.ones(3)*photons) )                
                self.polarTrianglesList.append(pTriangle)
    
        return np.array(self.trianglesList), np.array(self.polarTrianglesList)


    def generateFromRecArray(self, max_ast):
        self.trianglesList = []
        self.polarTrianglesList = []
        self.nfieldsSizes = []
        self.cumAstSizes = [0]
        for i in range(1,max_ast+1,1 ):
            number_of_asterisms = len(self.asterismsRecArray[i])
            print('number_of_asterisms:', number_of_asterisms)
            self.nfieldsSizes.append(number_of_asterisms)
            self.cumAstSizes.append(self.cumAstSizes[-1]+number_of_asterisms)
            for j in range(number_of_asterisms):
                xcoords = self.asterismsRecArray[i][j]['COORD'][0]
                ycoords = self.asterismsRecArray[i][j]['COORD'][1]
                flux = self.asterismsRecArray[i][j]['FLUXJ'] + self.asterismsRecArray[i][j]['FLUXH']
                triangle = np.vstack( [xcoords, ycoords, flux] )
                self.trianglesList.append(triangle)
                pTriangle = np.vstack( (cartesianToPolar(triangle[0:2,:]), flux) )                
                self.polarTrianglesList.append(pTriangle)                
        return np.array(self.trianglesList), np.array(self.polarTrianglesList)


    def getAsterismSource(self, field, ast, s):
        astIndexGlobal = self.cumAstSizes[field]+ast
        z = self.asterismsInputDataPolar[astIndexGlobal, 0, s]
        a = self.asterismsInputDataPolar[astIndexGlobal, 1, s]
        n = self.asterismsInputDataPolar[astIndexGlobal, 2, s]
        return np.array([z, a, n])


    def sourceIsPresent(self, source):
        # a source is a triple (array) of zenith, azimuth, photons
        for z, a, n, ii in zip(self.currentFieldSourcesData['Zenith'], 
                               self.currentFieldSourcesData['Azimuth'], 
                               self.currentFieldSourcesData['NumberPhotons'],
                               range(len(self.currentFieldSourcesData['Zenith'])) ):
            if np.allclose(np.array([z,a,n]), source):
                return ii
        return -1


    def getSourcesData(self, field):
        self.currentFieldSourcesData = {}
        self.currentFieldSourcesData['Zenith'] = [] 
        self.currentFieldSourcesData['Azimuth'] = []
        self.currentFieldSourcesData['NumberPhotons'] = []
        self.currentFieldAsterismsIndices = []        
        fieldsize = self.nfieldsSizes[field]
        for ast in range(fieldsize):
            astIndices = []
            for s in range(3):
                source = self.getAsterismSource(field, ast, s)
                s_index = self.sourceIsPresent(source)
                if s_index==-1:
                    self.currentFieldSourcesData['Zenith'].append(source[0])
                    self.currentFieldSourcesData['Azimuth'].append(source[1])
                    self.currentFieldSourcesData['NumberPhotons'].append(source[2])
                    s_index = len(self.currentFieldSourcesData['Zenith'])-1                    
                astIndices.append(s_index)                
            self.currentFieldAsterismsIndices.append(astIndices)


    def computeAsterisms(self, eeRadiusInMas):
        self.sr_Asterism = []
        self.fwhm_Asterism = []
        self.ee_Asterism = []
        self.cov_ellipses_Asterism = []
        for field in range(self.nfields):
            fieldsize = self.nfieldsSizes[field]
            print("*")            
            print("* Computing Field " + str(field) + " - number of asterisms :" + str(fieldsize))            
            print("*")            
            self.getSourcesData(field)
            base = self.cumAstSizes[field]
            for ast in range(fieldsize):
                print("*")
                print("* Computing Asterism " + str(ast))
                print("*")
                self.doOverallSimulation(ast)
                self.computeMetrics(eeRadiusInMas)
                self.sr_Asterism.append(self.sr) 
                self.fwhm_Asterism.append(self.fwhm) 
                self.ee_Asterism.append(self.ee) 
                self.cov_ellipses_Asterism.append(self.cov_ellipses)            
            covsarray = np.array(self.cov_ellipses_Asterism)[base:base+fieldsize, :,1]**2 + np.array(self.cov_ellipses_Asterism)[base:base+fieldsize, :,2]**2
            dd = np.average(np.abs(covsarray), axis=1)            
            self.twoPlots(dd, base)


    def plotAsterisms(self, dd, istart=None, istop=None, min_id=None):
        cmap = cm.get_cmap('rainbow')
        al = dd/600.0
        al = np.minimum(al, np.ones_like(dd))
        al = np.maximum(al, np.zeros_like(dd))
        al = 1.0 - al
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
        scales = 50.0 *  np.log(fluxes) 
        ntriangles = X.shape[0]
        plt.figure(figsize=(10, 10))
        plt.axvline(0, c='black', linewidth=0.5)
        plt.axhline(0, c='black', linewidth=0.5)        
        for i in range(istart, istop):
#            if i==min_id:
#                cc = [1.0, 0, 0]
#            else:
            cc = cmap(al[i-istart])
            aa = 1.0 # al[i-istart] * 0.9 + 0.1                
            coords = np.transpose(X[i, :2, :])
            xx = np.sum(X[i, 0, :])/3.0
            yy = np.sum(X[i, 1, :])/3.0
            px = X[i, 0, :]
            py = X[i, 1, :]
            plt.quiver([xx, xx, xx], [yy, yy, yy], px-xx, py-yy, color=cc, width=0.005, scale_units='xy', scale=1, alpha = aa )
#            t1 = plt.Polygon(coords), alpha = 0.1, color=( np.random.uniform(0, 1),
#                                                                             np.random.uniform(0, 1), 
#                                                                             np.random.uniform(0, 1)) )
#            plt.gca().add_patch(t1)

        coords = np.transpose(X[min_id, :2, :])
        t1 = plt.Polygon(coords, alpha = 0.3, color='r')
        
        plt.gca().add_patch(t1)
        plt.scatter(xcoords, ycoords, s=scales, c='yellow', edgecolors='r')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.gca().set_xlim([-60,60]) 
        plt.gca().set_ylim([-60,60]) 
        plt.show()

        
    def twoPlots(self, dd0, base):
        min_id = np.argmin(dd0)
        plt.plot(dd0)
        plt.show()
        self.plotAsterisms(dd0, base, base+dd0.shape[0], base+min_id)

