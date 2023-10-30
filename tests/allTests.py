from tiptop.tiptop import *
rc("text", usetex=False)

import unittest

def cpuArray(v):
    if isinstance(v,np.ndarray) or isinstance(v, list):
        return v
    else:
        return v.get()


class TestTiptop(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass
#        path = "data/ini/"
#        parametersFile = 'mavisParams'
#        fullPathFilename = path + parametersFile + '.ini'
#        windPsdFile = 'data/windpsd_mavis.fits'
#        TestMavisLO.mLO = MavisLO(path, parametersFile, verbose=True)
#        overallSimulation("perfTest", "SOUL", 'perfTest', 'testSOUL', doPlot=True, doConvolve=True)


class TestMavis(TestTiptop):

    def test_mavis(self):
        """
        Test
        """
        computed_result = overallSimulation("perfTest", "MAVIS", 'perfTest', 'testMAVIS', doPlot=False, doConvolve=True, returnRes=True)
        
        ii = 0 
#        for aa in computed_result:
#            with open('tests/mavisResult' + str(ii) +  '.npy', 'wb') as f:
#                np.save(f, cpuArray(aa))
#            ii += 1 
        
        stored_result0 = np.load('tests/mavisResult0.npy')
        stored_result1 = np.load('tests/mavisResult1.npy')

        self.assertTrue( np.testing.assert_allclose(cpuArray(computed_result[0]), stored_result0, rtol=1e-03, atol=1e-5)==None)
        self.assertTrue( np.testing.assert_allclose(cpuArray(computed_result[1]), stored_result1, rtol=1e-03, atol=1e-5)==None)


def suite():
    suite = unittest.TestSuite()
    suite.addTest(TestMavis('test_mavis'))
    return suite



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
