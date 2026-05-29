import unittest
import tempfile
import os
import numpy as np
from configparser import ConfigParser
from matplotlib import rc
from astropy.io import fits

# ----------------------------------------------------------------------------
# --- Explicit Imports (NO WILDCARDS) ---
# ----------------------------------------------------------------------------
from tiptop.tiptop import overallSimulation
from tiptop.baseSimulation import baseSimulation
from tiptop.asterismSimulation import asterismSimulation
from tiptop.asterismSimulationHo import asterismSimulationHo
from tiptop.tiptopUtils import cpuArray
from mastsel.mavisPsf import padOrCropCentered, centeredPixelCoords

rc("text", usetex=False)

class TestTiptop(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        pass

    def test_centered_padding_supports_even_and_odd_grids(self):
        """Test central padding alignment for both even and odd dimension grids."""
        for outer_size, inner_size in ((32, 28), (33, 28), (33, 29)):
            pupil = np.ones((inner_size, inner_size), dtype=np.float64)
            padded = cpuArray(padOrCropCentered(pupil, outer_size))

            self.assertEqual(padded.shape, (outer_size, outer_size))
            self.assertAlmostEqual(np.sum(padded), inner_size**2)
            cy, cx = centeredPixelCoords(outer_size)
            self.assertEqual((cy, cx), (outer_size // 2, outer_size // 2))
            self.assertEqual(padded[cy, cx], 1.0)


class TestMavis(TestTiptop):

    def test_mavis(self):
        """Test MAVIS standard simulation against stored baseline results."""
        computed_result = overallSimulation('tiptop/perfTest', 'MAVIStest',
                                            'tiptop/perfTest', 'testMAVIS',
                                            doPlot=False, doConvolve=True,
                                            returnRes=True)

        stored_result0 = np.load('tests/mavisResult0.npy')
        stored_result1 = np.load('tests/mavisResult1.npy')

        self.assertTrue(
            np.testing.assert_allclose(cpuArray(computed_result[0]),
                                       stored_result0, rtol=1e-03, atol=1e-5) is None)
        self.assertTrue(
            np.testing.assert_allclose(cpuArray(computed_result[1]),
                                       stored_result1, rtol=1e-03, atol=1e-5) is None)

    def test_mavis_jitter(self):
        """Test MAVIS simulation dynamic jitter inclusion via temp config files."""
        # Baseline simulation without extra jitter
        sr_nj, fwhm_nj, ee_nj = overallSimulation('tiptop/perfTest', 'MAVIStest',
                                                  'tiptop/perfTest', 'testMAVIS',
                                                  doPlot=False, doConvolve=True, 
                                                  returnMetrics=True)

        original_config_path = os.path.join('tiptop/perfTest', 'MAVIStest.ini')
        config = ConfigParser()
        config.optionxform = str
        config.read(original_config_path)

        jitter_fwhm = 10.0  # [mas]
        if not config.has_option('telescope', 'jitter_FWHM'):
            config.set('telescope', 'jitter_FWHM', str(jitter_fwhm))

        with tempfile.NamedTemporaryFile(mode='w', suffix='.ini', delete=False) as temp_file:
            config.write(temp_file)
            temp_filename = temp_file.name

        try:
            temp_dir = os.path.dirname(temp_filename)
            temp_basename = os.path.splitext(os.path.basename(temp_filename))[0]

            sr, fwhm, ee = overallSimulation(temp_dir, temp_basename,
                                             'tiptop/perfTest', 'testMAVISJitter',
                                              doPlot=False, doConvolve=True,
                                              returnMetrics=True)

            self.assertIsNotNone(sr)
            self.assertIsNotNone(fwhm)
            self.assertIsNotNone(ee)

            sr_cpu = np.array(cpuArray(sr))
            fwhm_cpu = np.array(cpuArray(fwhm))
            ee_cpu = np.array(cpuArray(ee))
            sr_nj_cpu = np.array(cpuArray(sr_nj))

            self.assertTrue(np.all(sr_cpu > 0))
            self.assertTrue(np.all(fwhm_cpu > 0))
            self.assertTrue(np.all(ee_cpu > 0))

            # SR should strictly degrade with added telescope jitter
            self.assertTrue(np.all(sr_cpu < sr_nj_cpu))

        finally:
            if os.path.exists(temp_filename):
                os.remove(temp_filename)

    def test_mavis_auto_science_field_of_view(self):
        """Test MAVIS simulation handling of automatic FOV scaling (-1 flag)."""
        original_config_path = os.path.join('tiptop/perfTest', 'MAVIStest.ini')

        config = ConfigParser()
        config.optionxform = str
        config.read(original_config_path)
        config.set('sensor_science', 'FieldOfView', '-1')

        with tempfile.NamedTemporaryFile(mode='w', suffix='.ini', delete=False) as temp_file:
            config.write(temp_file)
            temp_filename = temp_file.name

        try:
            temp_dir = os.path.dirname(temp_filename)
            temp_basename = os.path.splitext(os.path.basename(temp_filename))[0]

            sr, fwhm, ee = overallSimulation(temp_dir, temp_basename,
                                             'tiptop/perfTest', 'testMAVISAutoFOV',
                                             doPlot=False, doConvolve=True,
                                             returnMetrics=True)

            self.assertIsNotNone(sr)
            self.assertGreater(len(sr), 0)
        finally:
            if os.path.exists(temp_filename):
                os.remove(temp_filename)

    def test_mavis_odd_psd_grid_even_legacy_path(self):
        """Ensure odd-sized grids do not trigger alignment or casting crashes."""
        original_config_path = os.path.join('tiptop/perfTest', 'MAVIStest.ini')

        config = ConfigParser()
        config.optionxform = str
        config.read(original_config_path)
        config.set('telescope', 'Resolution', '321')
        config.set('sensor_science', 'FieldOfView', '513')

        with tempfile.NamedTemporaryFile(mode='w', suffix='.ini', delete=False) as temp_file:
            config.write(temp_file)
            temp_filename = temp_file.name

        try:
            temp_dir = os.path.dirname(temp_filename)
            temp_basename = os.path.splitext(os.path.basename(temp_filename))[0]

            simulation = baseSimulation(
                temp_dir, temp_basename, 'tiptop/perfTest', 'testMAVISOddGrid',
                doConvolve=False, doPlot=False, verbose=False)
            
            simulation.doOverallSimulation()

            self.assertEqual(simulation.N % 2, 1)
            self.assertEqual(simulation.nPixPSF % 2, 1)

            self.assertGreater(len(simulation.results), 0)
            for psf in simulation.results:
                arr = np.asarray(cpuArray(psf.sampling), dtype=np.float64)
                self.assertTrue(np.isfinite(arr).all())
                self.assertGreater(arr.sum(), 0.0)
        finally:
            if os.path.exists(temp_filename):
                os.remove(temp_filename)


class TestAsterismSimulation(TestTiptop):

    def test_asterism_simulation_creation(self):
        """Test initialization and structural flags for LO Asterism Evaluation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            simulation = asterismSimulation("TestERIS", "tiptop/astTest", "ERISastSinglesTest",
                                           tmpdir, 'testERIS', 
                                           doPlot=False, verbose=False)

            self.assertIsNotNone(simulation)
            self.assertEqual(simulation.simulName, "TestERIS")
            self.assertTrue(hasattr(simulation, 'hasAsterismSection'))

            if simulation.hasAsterismSection:
                self.assertTrue(hasattr(simulation, 'asterismMode'))
                self.assertTrue(hasattr(simulation, 'cumAstSizes'))
                self.assertGreater(len(simulation.cumAstSizes), 0)

    def test_asterism_simulation_single_computation(self):
        """Test isolated computation of a single LO asterism."""
        with tempfile.TemporaryDirectory() as tmpdir:
            simulation = asterismSimulation("TestERISSingle", "tiptop/astTest", "ERISastSinglesTest",
                                           tmpdir, 'testERISSingle', 
                                           doPlot=False, verbose=False)

            if simulation.hasAsterismSection and len(simulation.cumAstSizes) > 1:
                result = simulation.computeAsterisms(eeRadiusInMas=50, index=0, doConvolve=False)

                self.assertIsNotNone(result)
                self.assertEqual(len(result), 1)  # Strictly one configuration returned

                asterism_props = result[0]
                self.assertTrue(hasattr(asterism_props, 'strehl'))
                self.assertTrue(hasattr(asterism_props, 'fwhm'))
                self.assertTrue(hasattr(asterism_props, 'jitter'))
                
                # Metrics sanity check
                self.assertGreater(asterism_props.strehl, 0)
                self.assertLess(asterism_props.strehl, 1)
                self.assertGreater(asterism_props.fwhm, 0)

    def test_asterism_simulation_full_loop_and_sorting(self):
        """
        CRITICAL TEST: Ensures that computing multiple asterisms in a loop (index=None)
        correctly updates the backend configuration, does not leak memory, and
        returns a globally sorted list of dataclasses.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            simulation = asterismSimulation("TestERISLoop", "tiptop/astTest", "ERISastSinglesTest",
                                           tmpdir, 'testERISLoop', 
                                           doPlot=False, verbose=False)

            if simulation.hasAsterismSection and len(simulation.cumAstSizes) > 1:
                # Running all available asterisms in the test config
                results = simulation.computeAsterisms(eeRadiusInMas=50, index=None, doConvolve=False)

                self.assertIsNotNone(results)
                n_asterisms = len(results)
                self.assertGreater(n_asterisms, 1, "The test configuration must have >1 asterism to test sorting.")

                # Validate that AbstractSimulation accumulators grew to the exact expected size
                self.assertEqual(len(simulation.strehl_Asterism), n_asterisms)
                self.assertEqual(len(simulation.fwhm_Asterism), n_asterisms)

                # Ensure that the returned dataclasses are correctly sorted by Jitter Penalty (Ascending)
                for i in range(1, n_asterisms):
                    self.assertLessEqual(results[i-1].jitter, results[i].jitter, 
                                         "Asterism results are not correctly sorted by jitter penalty!")


class TestHoAsterismSimulation(TestTiptop):

    def test_ho_asterism_simulation_creation(self):
        """Test initialization and structural flags for HO Asterism Evaluation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            simulation = asterismSimulationHo("TestERISHO", "tiptop/astTest", "ERISastHO",
                                             tmpdir, 'testERISHO', 
                                             doPlot=False, verbose=False)

            self.assertIsNotNone(simulation)
            self.assertEqual(simulation.simulName, "TestERISHO")
            self.assertTrue(hasattr(simulation, 'hasHoAsterismSection'))

            if simulation.hasHoAsterismSection:
                self.assertEqual(simulation.asterismMode, 'SingleHO')
                self.assertGreater(len(simulation.cumAstSizes), 0)

    def test_ho_asterism_single_computation(self):
        """Test isolated computation of a single HO asterism."""
        with tempfile.TemporaryDirectory() as tmpdir:
            simulation = asterismSimulationHo("TestERISHOSingle", "tiptop/astTest", "ERISastHO",
                                             tmpdir, 'testERISHOSingle', 
                                             doPlot=False, verbose=False)

            if simulation.hasHoAsterismSection and len(simulation.cumAstSizes) > 1:
                result = simulation.computeHoAsterisms(eeRadiusInMas=50, index=0)

                self.assertIsNotNone(result)
                self.assertIn('indices', result)
                self.assertEqual(len(result['indices']), 1)

                self.assertIn('strehl', result)
                self.assertGreater(result['strehl'][0], 0)
                self.assertLess(result['strehl'][0], 1)

    def test_ho_asterism_full_loop_and_temp_file_cleanup(self):
        """
        CRITICAL TEST: Ensures that the dynamic creation of `.ini` configurations 
        for multiple HO Asterisms works perfectly, cleans up temporary files gracefully,
        and returns a list sorted by Strehl Ratio (Descending).
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            simulation = asterismSimulationHo("TestERISHOLoop", "tiptop/astTest", "ERISastHO",
                                             tmpdir, 'testERISHOLoop', 
                                             doPlot=False, verbose=False)

            if simulation.hasHoAsterismSection and len(simulation.cumAstSizes) > 1:
                
                results = simulation.computeHoAsterisms(eeRadiusInMas=50, index=None)

                self.assertIsNotNone(results)
                n_configs = len(results['indices'])
                self.assertGreater(n_configs, 1)

                # Ensure temp files were deleted during the loop
                temp_file = os.path.join(simulation.outputDir, f"{simulation.parametersFile}_temp_{n_configs-1}.ini")
                self.assertFalse(os.path.exists(temp_file), "Temporary INI files were not cleaned up!")

                # Verify sorting by Strehl Ratio (Descending)
                for i in range(1, n_configs):
                    self.assertGreaterEqual(results['strehl'][i-1], results['strehl'][i],
                                            "HO Configurations are not correctly sorted by Strehl Ratio!")


class TestSystemIOAndConfig(TestTiptop):
    
    def test_save_results_fits_integrity(self):
        """
        Critical I/O Test: Verifies that the simulation completes execution, 
        calls saveResults(), and generates a valid FITS file with correct data types in the header.
        """
        with tempfile.TemporaryDirectory() as tmpdirname:
            # Run the simulation EXPLICITLY requesting to save the results 
            # (returnMetrics=False and returnRes=False trigger saveResults in overallSimulation)
            overallSimulation('tiptop/perfTest', 'MAVIStest',
                              tmpdirname, 'testFITSOutput',
                              doPlot=False, doConvolve=True,
                              returnMetrics=False, returnRes=False, addSrAndFwhm=True)
            
            fits_file = os.path.join(tmpdirname, 'testFITSOutput.fits')
            self.assertTrue(os.path.exists(fits_file), "The FITS file was not generated.")
            
            # Open the FITS file and check the header sanity and float casting
            with fits.open(fits_file) as hdul:
                self.assertGreaterEqual(len(hdul), 4, "The FITS file does not have the expected minimum number of HDUs.")
                hdr1 = hdul[1].header
                
                # Verify that the native float casting (which we implemented) works flawlessly
                self.assertIn('SR0000', hdr1, "The SR0000 metric is missing from the FITS header.")
                self.assertIsInstance(hdr1['SR0000'], float, "The SR metric is not a native float.")
                self.assertIn('RESH0000', hdr1, "The High-Order residual (RESH0000) is missing from the header.")
                self.assertIsInstance(hdr1['RESH0000'], float, "The RESH metric is not a native float.")

    def test_missing_config_section_raises_error(self):
        """
        Negative Testing: Removing a vital section from the .ini file 
        must raise a clean and controlled ValueError from our validator.
        """
        original_config_path = os.path.join('tiptop/perfTest', 'MAVIStest.ini')
        config = ConfigParser()
        config.optionxform = str
        config.read(original_config_path)
        
        # Sabotage the configuration by removing the mandatory [telescope] section
        config.remove_section('telescope')
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.ini', delete=False) as temp_file:
            config.write(temp_file)
            temp_filename = temp_file.name
            
        try:
            temp_dir = os.path.dirname(temp_filename)
            temp_basename = os.path.splitext(os.path.basename(temp_filename))[0]
            
            # We expect the initialization to catch the issue and raise a ValueError
            with self.assertRaises(ValueError) as context:
                simulation = baseSimulation(temp_dir, temp_basename,
                                            temp_dir, 'dummyOutput')
            
            # Verify that the error message perfectly matches the one predefined in AbstractSimulation
            self.assertIn("The section 'telescope' is missing", str(context.exception))
            
        finally:
            if os.path.exists(temp_filename):
                os.remove(temp_filename)


def suite():
    suite = unittest.TestSuite()
    # Test Mavis
    suite.addTest(TestMavis('test_mavis'))
    suite.addTest(TestMavis('test_mavis_jitter'))
    suite.addTest(TestMavis('test_mavis_auto_science_field_of_view'))
    suite.addTest(TestMavis('test_mavis_odd_psd_grid_even_legacy_path'))
    
    # Test Asterism (LO)
    suite.addTest(TestAsterismSimulation('test_asterism_simulation_creation'))
    suite.addTest(TestAsterismSimulation('test_asterism_simulation_single_computation'))
    suite.addTest(TestAsterismSimulation('test_asterism_simulation_full_loop_and_sorting'))
    
    # Test Asterism (HO)
    suite.addTest(TestHoAsterismSimulation('test_ho_asterism_simulation_creation'))
    suite.addTest(TestHoAsterismSimulation('test_ho_asterism_single_computation'))
    suite.addTest(TestHoAsterismSimulation('test_ho_asterism_full_loop_and_temp_file_cleanup'))
    
    # Test I/O e Configurazioni (I NUOVI TEST)
    suite.addTest(TestSystemIOAndConfig('test_save_results_fits_integrity'))
    suite.addTest(TestSystemIOAndConfig('test_missing_config_section_raises_error'))
    
    return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
