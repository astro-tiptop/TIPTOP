from .baseSimulation import *
from dataclasses import dataclass
from typing import List
import itertools
from configparser import ConfigParser

@dataclass
class HoStar:
    zenith: float
    azimuth: float
    wavelength: float
    photons: float

@dataclass
class HoAsterismProperties:
    index: int
    ho_stars: List[HoStar]
    strehl_ratio: float
    fwhm: float
    encircled_energy: float
    ho_residual: float

def unrollHoAsterismData(all_combos, zenith, azimuth, wavelength, photons):
    """Unroll HO asterism data similar to the LO version"""
    asterism = np.array([np.take(zenith, all_combos),
                        np.take(azimuth, all_combos), 
                        np.take(wavelength, all_combos),
                        np.take(photons, all_combos)])
    return np.swapaxes(asterism, 0, 1)

class asterismSimulationHo(baseSimulation):

    def __init__(self, simulName, path, parametersFile, outputDir,
                 outputFile, doPlot=False, addSrAndFwhm=False, verbose=False,
                 getHoErrorBreakDown=False, progressStatus=False):

        # Initialize base simulation without LO part
        super().__init__(path, parametersFile, outputDir, outputFile, doConvolve=False,
                          doPlot=False, addSrAndFwhm=addSrAndFwhm,
                          verbose=verbose, getHoErrorBreakDown=getHoErrorBreakDown,
                          savePSDs=False)

        self.simulName = simulName
        self.doPlotAst = doPlot
        self.progressStatus = progressStatus
        self.firstConfigCall = True

        # HO asterism specific attributes
        self.asterismsInputDataHo = None
        self.hasHoAsterismSection = False
        self.nHoStars = 0

        # Results storage
        self.strehl_HoAsterism = []
        self.fwhm_HoAsterism = []
        self.ee_HoAsterism = []
        self.ho_res_HoAsterism = []

        if 'HO_ASTERISM_SELECTION' in self.my_data_map.keys():
            self.hasHoAsterismSection = True
            self.asterismMode = self.my_data_map['HO_ASTERISM_SELECTION']['mode']

            # Read HO stars configuration
            listZ = self.my_data_map['HO_ASTERISM_SELECTION']['Zenith']
            listA = self.my_data_map['HO_ASTERISM_SELECTION']['Azimuth'] 
            listW = self.my_data_map['HO_ASTERISM_SELECTION']['Wavelength']
            listP = self.my_data_map['HO_ASTERISM_SELECTION']['NumberPhotons']

            self.nHoStars = len(listZ)
            self.cumAstSizes = [0]
            self.nfields = 1

            if self.asterismMode == 'SingleHO':
                # Each HO star is evaluated individually
                all_combos = list(itertools.combinations(list(range(self.nHoStars)), 1))
                self.nfieldsSizes = [len(all_combos)]
                self.cumAstSizes.append(self.nfieldsSizes[0])

                # Convert to numpy arrays
                zenith = np.array(listZ, dtype=np.float64)
                azimuth = np.array(listA, dtype=np.float64)
                wavelength = np.array(listW, dtype=np.float64)
                photons = np.array(listP, dtype=np.float64)

                self.asterismsInputDataHo = unrollHoAsterismData(all_combos, zenith, azimuth, wavelength, photons)
                self.allHoAsterismsIndices = np.asarray(all_combos)

            else:
                raise ValueError(f"Unknown HO asterism mode: {self.asterismMode}")

            if self.verbose:
                print(f'HO Asterism mode: {self.asterismMode}')
                print(f'Number of HO configurations: {self.cumAstSizes[-1]}')

    def configHO(self, hoAsterismIndex):
        """Configure HO sources for a specific asterism"""
        if hoAsterismIndex is None:
            return

        # Get the HO star indices for this asterism
        ho_star_indices = self.allHoAsterismsIndices[hoAsterismIndex]

        # Extract HO configuration for this asterism
        ho_asterism_data = self.asterismsInputDataHo[hoAsterismIndex]

        # Update sources_HO in the configuration
        if len(ho_asterism_data[0]) == 1:
            self.my_data_map['sources_HO']['Zenith'] = [float(ho_asterism_data[0][0])]
            self.my_data_map['sources_HO']['Azimuth'] = [float(ho_asterism_data[1][0])]
            self.my_data_map['sources_HO']['Wavelength'] = [float(ho_asterism_data[2][0])]
            self.my_data_map['sensor_HO']['NumberPhotons'] = [float(ho_asterism_data[3][0])]
        else:
            self.my_data_map['sources_HO']['Zenith'] = ho_asterism_data[0].tolist()
            self.my_data_map['sources_HO']['Azimuth'] = ho_asterism_data[1].tolist()
            self.my_data_map['sources_HO']['Wavelength'] = ho_asterism_data[2].tolist()
            self.my_data_map['sensor_HO']['NumberPhotons'] = ho_asterism_data[3].tolist()

        # Write a temporary configuration file for this asterism
        temp_filename = f"{self.parametersFile}_temp_{hoAsterismIndex}"
        temp_fullpath = os.path.join(self.outputDir, temp_filename + '.ini')

        try:
            # Write the temporary file
            config = ConfigParser()
            config.optionxform = str
            for section_name, section_data in self.my_data_map.items():
                if section_name == 'HO_ASTERISM_SELECTION':
                    continue
                config.add_section(section_name)
                for key, value in section_data.items():
                    if isinstance(value, str):
                        config.set(section_name, key, f"'{value}'")
                    else:
                        config.set(section_name, key, str(value))

            with open(temp_fullpath, 'w') as f:
                config.write(f)

            # Update the path for fourierModel
            self.temp_parametersFile = temp_filename
            self.temp_path = self.outputDir

            if self.verbose:
                print(f'Configured HO asterism {hoAsterismIndex} with {len(ho_star_indices)} stars')
                print(f'Temporary config file: {temp_fullpath}')
                for i, idx in enumerate(ho_star_indices):
                    print(f'  HO Star {i}: Zenith={ho_asterism_data[0][i]:.2f}, Azimuth={ho_asterism_data[1][i]:.2f}, '
                        f'Wavelength={ho_asterism_data[2][i]*1e9:.0f}nm, Photons={ho_asterism_data[3][i]:.1f}')

        except Exception as e:
            print(f"Error creating temporary configuration file: {e}")
            # Remove temp attributes if file creation failed
            if hasattr(self, 'temp_parametersFile'):
                delattr(self, 'temp_parametersFile')
            if hasattr(self, 'temp_path'):
                delattr(self, 'temp_path')
            raise


    def computeHoAsterisms(self, eeRadiusInMas=50, index=None):
        """Compute HO asterisms performance"""

        if index is None:
            singleAsterism = False
            nConfigs = self.nfieldsSizes[0]
            configs_to_process = range(nConfigs)
        else:
            singleAsterism = True
            configs_to_process = [index]

        self.eeRadiusInMas = eeRadiusInMas
        self.strehl_HoAsterism = []
        self.fwhm_HoAsterism = []
        self.ee_HoAsterism = []
        self.ho_res_HoAsterism = []

        for config_idx in configs_to_process:
            if self.progressStatus:
                print(f'Processing HO configuration {config_idx+1}/{len(configs_to_process)}')

            # Configure HO sources for this configuration
            self.configHO(config_idx)

            # Reset first call flag to force recalculation
            self.firstSimCall = True

            # Store original values
            original_path = self.path
            original_parametersFile = self.parametersFile
            original_fullPathFilename = self.fullPathFilename

            if hasattr(self, 'temp_path'):
                print(f'Using temporary path: {self.temp_path}, file: {self.temp_parametersFile}')
                try:
                    # Load configuration file
                    self.loadConfigurationFile(path=self.temp_path, parametersFile=self.temp_parametersFile)
                except Exception as e:
                    print(f"Error loading temporary configuration: {e}")
                    # Restore original values and continue with next config
                    self.path = original_path
                    self.parametersFile = original_parametersFile
                    self.fullPathFilename = original_fullPathFilename
                    continue
            else:
                print('Temporary path or filename not set. Using original configuration.')

            try:
                # Run the simulation for this HO configuration
                self.doOverallSimulation()

                # Compute and store metrics
                self.computeMetrics()

                # Store results
                self.strehl_HoAsterism.append(np.array([cpuArray(x) for x in self.sr]))
                self.fwhm_HoAsterism.append(self.fwhm)
                self.ee_HoAsterism.append(self.ee)
                self.ho_res_HoAsterism.append(np.array([cpuArray(x) for x in self.HO_res]))

                if self.verbose:
                    print(f'Config {config_idx}: SR={self.sr[0]:.4f}, FWHM={self.fwhm[0]:.2f}mas, '
                        f'EE={self.ee[0]:.4f}, HO_res={self.HO_res[0]:.1f}nm')
   
            except Exception as e:
                print(f"Error in simulation for config {config_idx}: {e}")
                # Store NaN values for failed simulations
                self.strehl_HoAsterism.append(np.array([np.nan]))
                self.fwhm_HoAsterism.append([np.nan])
                self.ee_HoAsterism.append([np.nan])
                self.ho_res_HoAsterism.append(np.array([np.nan]))

            finally:
                # Restore original path and filename
                self.path = original_path
                self.parametersFile = original_parametersFile
                self.fullPathFilename = original_fullPathFilename

                # Clean up temporary file
                if hasattr(self, 'temp_parametersFile'):
                    temp_file = os.path.join(self.temp_path, self.temp_parametersFile + '.ini')
                    if os.path.exists(temp_file):
                        try:
                            os.remove(temp_file)
                        except Exception as e:
                            if self.verbose:
                                print(f"Warning: Could not remove temporary file {temp_file}: {e}")

            if self.verbose:
                print(f'Config {config_idx}: SR={self.sr[0]:.4f}, FWHM={self.fwhm[0]:.2f}mas, '
                      f'EE={self.ee[0]:.4f}, HO_res={self.HO_res[0]:.1f}nm')

        if not singleAsterism:
            # Save results
            np.save(os.path.join(self.outputDir, self.simulName+'_ho_sr.npy'), np.array(self.strehl_HoAsterism))
            np.save(os.path.join(self.outputDir, self.simulName+'_ho_fw.npy'), np.array(self.fwhm_HoAsterism))
            np.save(os.path.join(self.outputDir, self.simulName+'_ho_ee.npy'), np.array(self.ee_HoAsterism))
            np.save(os.path.join(self.outputDir, self.simulName+'_ho_res.npy'), np.array(self.ho_res_HoAsterism))

        if singleAsterism:
            # Return single result
            ho_stars_data = self.asterismsInputDataHo[index]
            ho_stars = []
            for i in range(len(ho_stars_data[0])):
                ho_stars.append(HoStar(ho_stars_data[0][i], ho_stars_data[1][i], 
                                      ho_stars_data[2][i], ho_stars_data[3][i]))

            return [HoAsterismProperties(index, ho_stars, self.sr[0], self.fwhm[0], 
                                       self.ee[0], self.HO_res[0])]
        else:
            # Return all results
            results = []
            sorted_indices = np.argsort([sr[0] for sr in self.strehl_HoAsterism])[::-1]  # Sort by SR descending

            for i, idx in enumerate(sorted_indices):
                ho_stars_data = self.asterismsInputDataHo[idx]
                ho_stars = []
                for j in range(len(ho_stars_data[0])):
                    ho_stars.append(HoStar(ho_stars_data[0][j], ho_stars_data[1][j],
                                          ho_stars_data[2][j], ho_stars_data[3][j]))

                results.append(HoAsterismProperties(idx, ho_stars, 
                                                  self.strehl_HoAsterism[idx][0],
                                                  self.fwhm_HoAsterism[idx][0],
                                                  self.ee_HoAsterism[idx][0],
                                                  self.ho_res_HoAsterism[idx][0]))

            return results

    def reloadHoResults(self):
        """Reload previously computed results"""
        self.strehl_HoAsterism = np.load(os.path.join(self.outputDir, self.simulName+'_ho_sr.npy'))
        self.fwhm_HoAsterism = np.load(os.path.join(self.outputDir, self.simulName+'_ho_fw.npy'))
        self.ee_HoAsterism = np.load(os.path.join(self.outputDir, self.simulName+'_ho_ee.npy'))
        self.ho_res_HoAsterism = np.load(os.path.join(self.outputDir, self.simulName+'_ho_res.npy'))

    def plotHoResults(self):
        """Plot HO asterism results"""
        if len(self.strehl_HoAsterism) == 0:
            print("No results to plot. Run computeHoAsterisms first.")
            return

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # Strehl Ratio
        sr_values = [sr[0] for sr in self.strehl_HoAsterism]
        axes[0,0].bar(range(len(sr_values)), sr_values)
        axes[0,0].set_title('Strehl Ratio')
        axes[0,0].set_xlabel('Configuration Index')
        axes[0,0].set_ylabel('Strehl Ratio')

        # FWHM
        fwhm_values = [fwhm[0] for fwhm in self.fwhm_HoAsterism]
        axes[0,1].bar(range(len(fwhm_values)), fwhm_values)
        axes[0,1].set_title('FWHM [mas]')
        axes[0,1].set_xlabel('Configuration Index')
        axes[0,1].set_ylabel('FWHM [mas]')

        # Encircled Energy
        ee_values = [ee[0] for ee in self.ee_HoAsterism]
        axes[1,0].bar(range(len(ee_values)), ee_values)
        axes[1,0].set_title(f'Encircled Energy @ {self.eeRadiusInMas} mas')
        axes[1,0].set_xlabel('Configuration Index')
        axes[1,0].set_ylabel('Encircled Energy')

        # HO Residual
        ho_res_values = [ho_res[0] for ho_res in self.ho_res_HoAsterism]
        axes[1,1].bar(range(len(ho_res_values)), ho_res_values)
        axes[1,1].set_title('HO Residual [nm RMS]')
        axes[1,1].set_xlabel('Configuration Index')
        axes[1,1].set_ylabel('HO Residual [nm RMS]')

        plt.tight_layout()
        plt.show()

        # Print summary
        best_sr_idx = np.argmax(sr_values)
        print(f"\nBest configuration (highest SR): {best_sr_idx}")
        print(f"  Strehl Ratio: {sr_values[best_sr_idx]:.4f}")
        print(f"  FWHM: {fwhm_values[best_sr_idx]:.2f} mas")
        print(f"  EE: {ee_values[best_sr_idx]:.4f}")
        print(f"  HO Residual: {ho_res_values[best_sr_idx]:.1f} nm RMS")