from tiptop.tiptop import *
rc("text", usetex=False)

import unittest
import tempfile
import os
from configparser import ConfigParser


def cpuArray(v):
    if isinstance(v,np.ndarray) or isinstance(v, list):
        return v
    else:
        return v.get()


class TestTiptop(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass
#        path = "tiptop/data/ini/"
#        parametersFile = 'mavisParams'
#        fullPathFilename = path + parametersFile + '.ini'
#        windPsdFile = 'tiptop/data/windpsd_mavis.fits'
#        TestMavisLO.mLO = MavisLO(path, parametersFile, verbose=True)
#        overallSimulation("perfTest", "SOUL", 'perfTest', 'testSOUL', doPlot=True, doConvolve=True)


class TestMavis(TestTiptop):

    def test_mavis(self):
        """
        Test MAVIS simulation against stored results
        """
        computed_result = overallSimulation('tiptop/perfTest', 'MAVIStest', 'tiptop/perfTest', 'testMAVIS',
                                            doPlot=False, doConvolve=True, returnRes=True)

        # This can be used to save / update the results if needed
        save_results = False
        if save_results:
            ii = 0
            for aa in computed_result:
                with open('tests/mavisResult' + str(ii) +  '.npy', 'wb') as f:
                    np.save(f, cpuArray(aa))
                ii += 1

        stored_result0 = np.load('tests/mavisResult0.npy')
        stored_result1 = np.load('tests/mavisResult1.npy')

        self.assertTrue( np.testing.assert_allclose(cpuArray(computed_result[0]), stored_result0, rtol=1e-03, atol=1e-5)==None)
        self.assertTrue( np.testing.assert_allclose(cpuArray(computed_result[1]), stored_result1, rtol=1e-03, atol=1e-5)==None)

    def test_mavis_jitter(self):
        """
        Test MAVIS simulation with jitter_FWHM
        """

        # Run simulation with temporary file
        sr_nj, fwhm_nj, ee_nj = overallSimulation('tiptop/perfTest', 'MAVIStest', 'tiptop/perfTest', 'testMAVIS',
                                            doPlot=False, doConvolve=True, returnMetrics=True)

        # Load the base configuration
        original_config_path = os.path.join('tiptop/perfTest', 'MAVIStest.ini')

        # Read the original configuration
        config = ConfigParser()
        config.optionxform = str
        config.read(original_config_path)

        jitter_fwhm = 10.0  # Example jitter FWHM value in mas

        # Add jitter_FWHM to the telescope section
        if not config.has_option('telescope', 'jitter_FWHM'):
            config.set('telescope', 'jitter_FWHM', str(jitter_fwhm))

        # Create a temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.ini', delete=False) as temp_file:
            config.write(temp_file)
            temp_filename = temp_file.name

        try:
            # Extract filename without extension and directory
            temp_dir = os.path.dirname(temp_filename)
            temp_basename = os.path.splitext(os.path.basename(temp_filename))[0]

            # Run simulation with temporary file
            sr, fwhm, ee = overallSimulation(temp_dir, temp_basename, 'tiptop/perfTest', 'testMAVISJitter',
                                              doPlot=False, doConvolve=True, returnMetrics=True)

            # Verify that the result is valid
            self.assertIsNotNone(sr)
            self.assertIsNotNone(fwhm)
            self.assertIsNotNone(ee)

            # Verify that values are reasonable (non-zero)
            self.assertGreater(len(sr), 0)  # sr should have elements
            self.assertGreater(len(fwhm), 0)  # fwhm should have elements
            self.assertGreater(len(ee), 0)  # ee should have elements

            sr_cpu = np.array(cpuArray(sr))
            fwhm_cpu = np.array(cpuArray(fwhm))
            ee_cpu = np.array(cpuArray(ee))

            sr_nj_cpu = np.array(cpuArray(sr_nj))
            fwhm_nj_cpu = np.array(cpuArray(fwhm_nj))

            # Verify that values are numeric and positive
            self.assertTrue(np.all(sr_cpu > 0))
            self.assertTrue(np.all(fwhm_cpu > 0))
            self.assertTrue(np.all(ee_cpu > 0))

            # SR should be lower with jitter
            self.assertTrue(np.all(sr_cpu < sr_nj_cpu))

        finally:
            # Cleanup: remove temporary file
            if os.path.exists(temp_filename):
                os.remove(temp_filename)


class TestAsterismSimulation(TestTiptop):

    def test_asterism_simulation_creation(self):
        """
        Test that the asterismSimulation class initializes correctly
        """
        # Test with ERIS standard configuration
        simulation = asterismSimulation("TestERIS", "tiptop/astTest", "ERISastSinglesTest",
                                       'tiptop/astTest', 'testERIS', 
                                       doPlot=False, verbose=False)

        # Check that the simulation was created correctly
        self.assertIsNotNone(simulation)
        self.assertEqual(simulation.simulName, "TestERIS")
        self.assertTrue(hasattr(simulation, 'hasAsterismSection'))

        # If there is an asterism section, check the basic parameters
        if simulation.hasAsterismSection:
            self.assertTrue(hasattr(simulation, 'asterismMode'))
            self.assertTrue(hasattr(simulation, 'cumAstSizes'))
            self.assertGreater(len(simulation.cumAstSizes), 0)

    def test_asterism_simulation_single_computation(self):
        """
        Test calculation of a single asterism
        """
        # Use ERIS configuration that should be available
        simulation = asterismSimulation("TestERISSingle", "tiptop/astTest", "ERISastSinglesTest",
                                       'tiptop/astTest', 'testERISSingle', 
                                       doPlot=False, verbose=False)

        if simulation.hasAsterismSection and len(simulation.cumAstSizes) > 1:
            # Test the calculation of a single asterism (index 0)
            result = simulation.computeAsterisms(eeRadiusInMas=50, index=0, doConvolve=False)

            # Check that the result is valid
            self.assertIsNotNone(result)
            self.assertEqual(len(result), 1)  # One asterism

            # Check the properties of the result
            asterism_props = result[0]
            self.assertTrue(hasattr(asterism_props, 'strehl'))
            self.assertTrue(hasattr(asterism_props, 'fwhm'))
            self.assertTrue(hasattr(asterism_props, 'jitter'))

            # Check that the values are reasonable
            self.assertGreater(asterism_props.strehl, 0)
            self.assertLess(asterism_props.strehl, 1)
            self.assertGreater(asterism_props.fwhm, 0)
            self.assertGreater(asterism_props.jitter, 0)


class TestHoAsterismSimulation(TestTiptop):

    def test_ho_asterism_simulation_creation(self):
        """
        Test that the asterismSimulationHo class initializes correctly
        """
        # Test with HO configuration
        simulation = asterismSimulationHo("TestERISHO", "tiptop/astTest", "ERISastHO",
                                         'tiptop/astTest', 'testERISHO', 
                                         doPlot=False, verbose=False)

        # Check that the simulation was created correctly
        self.assertIsNotNone(simulation)
        self.assertEqual(simulation.simulName, "TestERISHO")
        self.assertTrue(hasattr(simulation, 'hasHoAsterismSection'))

        # If there is a HO asterism section, check the basic parameters
        if simulation.hasHoAsterismSection:
            self.assertTrue(hasattr(simulation, 'asterismMode'))
            self.assertTrue(hasattr(simulation, 'cumAstSizes'))
            self.assertTrue(hasattr(simulation, 'asterismsInputDataHo'))
            self.assertGreater(len(simulation.cumAstSizes), 0)
            self.assertEqual(simulation.asterismMode, 'SingleHO')

    def test_ho_asterism_single_computation(self):
        """
        Test the calculation of a single HO configuration
        """
        simulation = asterismSimulationHo("TestERISHOSingle", "tiptop/astTest", "ERISastHO",
                                         'tiptop/astTest', 'testERISHOSingle', 
                                         doPlot=False, verbose=False)

        if simulation.hasHoAsterismSection and len(simulation.cumAstSizes) > 1:
            # Test the calculation of a single HO configuration (index 0)
            result = simulation.computeHoAsterisms(eeRadiusInMas=50, index=0)

            # Check that the result is valid
            self.assertIsNotNone(result)
            self.assertEqual(len(result), 1)  # One configuration

            # Check the properties of the result
            ho_props = result[0]
            self.assertTrue(hasattr(ho_props, 'strehl_ratio'))
            self.assertTrue(hasattr(ho_props, 'fwhm'))
            self.assertTrue(hasattr(ho_props, 'ho_residual'))
            self.assertTrue(hasattr(ho_props, 'ho_stars'))

            # Check that the values are reasonable
            self.assertGreater(ho_props.strehl_ratio, 0)
            self.assertLess(ho_props.strehl_ratio, 1)
            self.assertGreater(ho_props.fwhm, 0)
            self.assertGreater(ho_props.ho_residual, 0)
            self.assertGreater(len(ho_props.ho_stars), 0)

    def test_ho_asterism_config_temp_file(self):
        """
        Test that the creation of temporary files works correctly
        """
        simulation = asterismSimulationHo("TestERISHOTemp", "tiptop/astTest", "ERISastHO",
                                         'tiptop/astTest', 'testERISHOTemp', 
                                         doPlot=False, verbose=False)

        if simulation.hasHoAsterismSection:
            # Test the HO configuration (which creates temporary files)
            simulation.configHO(0)

            # Check that the temporary attributes have been created
            self.assertTrue(hasattr(simulation, 'temp_parametersFile'))
            self.assertTrue(hasattr(simulation, 'temp_path'))

            # Check that the temporary file exists
            temp_file = os.path.join(simulation.temp_path, simulation.temp_parametersFile + '.ini')
            self.assertTrue(os.path.exists(temp_file))

            # Cleanup
            if os.path.exists(temp_file):
                os.remove(temp_file)


def suite():
    suite = unittest.TestSuite()
    suite.addTest(TestMavis('test_mavis'))
    suite.addTest(TestAsterismSimulation('test_asterism_simulation_creation'))
    suite.addTest(TestAsterismSimulation('test_asterism_simulation_single_computation'))
    suite.addTest(TestHoAsterismSimulation('test_ho_asterism_simulation_creation'))
    suite.addTest(TestHoAsterismSimulation('test_ho_asterism_single_computation'))
    suite.addTest(TestHoAsterismSimulation('test_ho_asterism_config_temp_file'))
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())