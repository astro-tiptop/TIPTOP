import os
import itertools
import numpy as np
import matplotlib.pyplot as plt
from configparser import ConfigParser

# TIPTOP explicit base imports
from .baseSimulation import baseSimulation
from .tiptopUtils import cpuArray


def unrollHoAsterismData(all_combos, zenith, azimuth, wavelength, photons):
    """Unroll High-Order asterism parameters into indexed configurations."""
    asterism = np.array([np.take(zenith, all_combos),
                        np.take(azimuth, all_combos), 
                        np.take(wavelength, all_combos),
                        np.take(photons, all_combos)])
    return np.swapaxes(asterism, 0, 1)


class asterismSimulationHo(baseSimulation):
    """
    Evaluates individual or segmented High-Order configurations using P3 backend logic.
    Inherits lifecycle tracking from baseSimulation.
    """

    def __init__(self, simulName, path, parametersFile, outputDir,
                 outputFile, doPlot=False, addSrAndFwhm=False, verbose=False,
                 getHoErrorBreakDown=False, progressStatus=False):

        # Base structure is initialized with zero low order additions
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
        """Generates dynamic local ini configurations for segmented high-order sources."""
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
                for i, idx in enumerate(ho_star_indices):
                    print(f'  HO Star {i}: Zenith={ho_asterism_data[0][i]:.2f}, Azimuth={ho_asterism_data[1][i]:.2f}, Wavelength={ho_asterism_data[2][i]*1e9:.0f}nm')

        except Exception as e:
            print(f"Error creating temporary configuration file: {e}")
            if hasattr(self, 'temp_parametersFile'): delattr(self, 'temp_parametersFile')
            if hasattr(self, 'temp_path'): delattr(self, 'temp_path')
            raise


    def computeHoAsterisms(self, eeRadiusInMas=50, index=None):
        """Loops across the set configurations, resetting parameters safely between evaluations."""
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
                if self.verbose:
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
                self.doOverallSimulation(astIndex=None)
                self.computeMetrics()

                # --- SANITIZE AND HOMOGENIZE METRICS (LGS) ---
                # Extract arrays, flatten them, and cast each element to a native float
                clean_sr = [float(x) for x in np.atleast_1d(np.squeeze(cpuArray(self.sr)))]
                clean_fwhm = [float(x) for x in np.atleast_1d(np.squeeze(cpuArray(self.fwhm)))]
                clean_ee = [float(x) for x in np.atleast_1d(np.squeeze(cpuArray(self.ee)))]
                clean_ho_res = [float(x) for x in np.atleast_1d(np.squeeze(cpuArray(self.HO_res)))]

                # In HO mode, global metrics (SR, EE) might be returned as single scalars. 
                # We duplicate them to match the number of pointings (FWHM) to ensure 
                # the exact same dimensionality as the NGS case.
                if len(clean_sr) == 1 and len(clean_fwhm) > 1:
                    clean_sr = clean_sr * len(clean_fwhm)
                if len(clean_ee) == 1 and len(clean_fwhm) > 1:
                    clean_ee = clean_ee * len(clean_fwhm)

                self.strehl_HoAsterism.append(clean_sr)
                self.fwhm_HoAsterism.append(clean_fwhm)
                self.ee_HoAsterism.append(clean_ee)
                self.ho_res_HoAsterism.append(clean_ho_res)

                if self.verbose:
                    print(f'Config {config_idx}: SR={self.sr[0]:.4f}, FWHM={self.fwhm[0]:.2f}mas')
   
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
                            if self.verbose: print(f"Warning: Could not remove temporary file: {e}")

        if not singleAsterism:
            # Save results
            np.save(os.path.join(self.outputDir, self.simulName+'_ho_sr.npy'), np.array(self.strehl_HoAsterism))
            np.save(os.path.join(self.outputDir, self.simulName+'_ho_fw.npy'), np.array(self.fwhm_HoAsterism))
            np.save(os.path.join(self.outputDir, self.simulName+'_ho_ee.npy'), np.array(self.ee_HoAsterism))
            np.save(os.path.join(self.outputDir, self.simulName+'_ho_res.npy'), np.array(self.ho_res_HoAsterism))

        if singleAsterism:
            # Return single asterism data packed in numpy arrays
            ho_stars_data = self.asterismsInputDataHo[index]
            return {
                'indices': np.array([index]),
                'zenith': np.array([ho_stars_data[0]]),
                'azimuth': np.array([ho_stars_data[1]]),
                'wavelength': np.array([ho_stars_data[2]]),
                'photons': np.array([ho_stars_data[3]]),
                'strehl': np.array([cpuArray(self.sr[0])]),
                'fwhm': np.array([cpuArray(self.fwhm[0])]),
                'ee': np.array([cpuArray(self.ee[0])]),
                'ho_res': np.array([cpuArray(self.HO_res[0])])
            }
        else:
            # Extract Strehl Ratios and sort indices in descending Strehl order
            strehls = np.array([sr[0] for sr in self.strehl_HoAsterism])
            sorted_indices = np.argsort(strehls)[::-1]

            # Apply the sorting to the original data matrix in one step
            sorted_ho_data = self.asterismsInputDataHo[sorted_indices]

            # Constructs the output dictionary fully vectorized
            return {
                'indices': sorted_indices,
                'zenith': sorted_ho_data[:, 0, :],       # Shape: (Num_Ast, Stars_Per_Ast)
                'azimuth': sorted_ho_data[:, 1, :],
                'wavelength': sorted_ho_data[:, 2, :],
                'photons': sorted_ho_data[:, 3, :],
                'strehl': strehls[sorted_indices],
                'fwhm': np.array([self.fwhm_HoAsterism[idx][0] for idx in sorted_indices]),
                'ee': np.array([self.ee_HoAsterism[idx][0] for idx in sorted_indices]),
                'ho_res': np.array([self.ho_res_HoAsterism[idx][0] for idx in sorted_indices])
            }


    def reloadHoResults(self):
        """Reload previously computed results"""
        self.strehl_HoAsterism = np.load(os.path.join(self.outputDir, self.simulName+'_ho_sr.npy')).tolist()
        self.fwhm_HoAsterism = np.load(os.path.join(self.outputDir, self.simulName+'_ho_fw.npy')).tolist()
        self.ee_HoAsterism = np.load(os.path.join(self.outputDir, self.simulName+'_ho_ee.npy')).tolist()
        self.ho_res_HoAsterism = np.load(os.path.join(self.outputDir, self.simulName+'_ho_res.npy')).tolist()


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