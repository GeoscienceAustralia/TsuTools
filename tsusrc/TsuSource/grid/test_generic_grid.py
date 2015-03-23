"""Test suite for generic_grid.py
"""

import unittest
import os
import numpy
import generic_grid

# Check that GDAL interface is installed
try:
    import gdal
    import osr
    gdal_installed = True
except ImportError:
    gdal_installed = False

def test_function(x,y):
    return x*x + 2.5*y

class Test_generic_grid(unittest.TestCase):
    """
    Tests that grids are created correctly and for evaluation
    of functions over a grid
    """

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_create_grid(self):
        test_data={'xmin': 0, 'xmax': 5,
                   'ymin': -5, 'ymax':3, 
                   'dx': 1, 'dy': 1,
                   'xx_test': numpy.array([[0,1,2,3,4],[0,1,2,3,4],[0,1,2,3,4],
                                         [0,1,2,3,4],[0,1,2,3,4],[0,1,2,3,4],
                                         [0,1,2,3,4],[0,1,2,3,4]]),
                   'yy_test': numpy.array([[-5,-5,-5,-5,-5],[-4,-4,-4,-4,-4],
                                         [-3,-3,-3,-3,-3],[-2,-2,-2,-2,-2],
                                         [-1,-1,-1,-1,-1],[0,0,0,0,0],
                                         [1,1,1,1,1],[2,2,2,2,2]])}

        x, y, xx, yy = generic_grid.create_grid(test_data['xmin'],
                                                test_data['xmax'],
                                                test_data['ymin'],
                                                test_data['ymax'],
                                                test_data['dx'],
                                                test_data['dy'])

        msg = 'Grided x values not as expected in create_grid()'
        assert numpy.allclose(test_data['xx_test'], xx, rtol=1.e-9), msg
        msg = 'Grided y values not as expected in create_grid()'
        assert numpy.allclose(test_data['yy_test'], yy, rtol=1.e-9), msg
                            
    def test_evaluate_function_over_grid(self):

        test_data={'xmin': 0, 'xmax': 5,
                   'ymin': -5, 'ymax':3, 
                   'dx': 1, 'dy': 1,
                   'xx': numpy.array([[0,1,2,3,4],[0,1,2,3,4],[0,1,2,3,4],
                                         [0,1,2,3,4],[0,1,2,3,4],[0,1,2,3,4],
                                         [0,1,2,3,4],[0,1,2,3,4]]),
                   'yy': numpy.array([[-5,-5,-5,-5,-5],[-4,-4,-4,-4,-4],
                                         [-3,-3,-3,-3,-3],[-2,-2,-2,-2,-2],
                                         [-1,-1,-1,-1,-1],[0,0,0,0,0],
                                         [1,1,1,1,1],[2,2,2,2,2]]),
                   'z': numpy.array([[-12.5, -11.5,  -8.5,  -3.5,   3.5],
                                     [-10. ,  -9. ,  -6. ,  -1. ,   6. ],
                                     [ -7.5,  -6.5,  -3.5,   1.5,   8.5],
                                     [ -5. ,  -4. ,  -1. ,   4. ,  11. ],
                                     [ -2.5,  -1.5,   1.5,   6.5,  13.5],
                                     [  0. ,   1. ,   4. ,   9. ,  16. ],
                                     [  2.5,   3.5,   6.5,  11.5,  18.5],
                                     [  5. ,   6. ,   9. ,  14. ,  21. ]])}


        # Test 2 different ways of calling the function
        x, y, xx, yy, z = generic_grid.evaluate_function_over_grid(test_function,
                                                                  gridx=test_data['xx'],
                                                                  gridy=test_data['yy'])
        msg = 'Calculated values different to expect values in' \
            'evaluate_function_over_grid'
        assert numpy.allclose(test_data['z'], z, rtol=1.e-9), msg
        x, y, xx, yy, z = generic_grid.evaluate_function_over_grid(test_function,
                                                             xmin=test_data['xmin'],
                                                             xmax=test_data['xmax'],
                                                             ymin=test_data['ymin'],
                                                             ymax=test_data['ymax'],
                                                             dx=test_data['dx'],
                                                             dy=test_data['dy'])
        msg = 'Calculated values different to expect values in' \
            'evaluate_function_over_grid'
        assert numpy.allclose(test_data['z'], z, rtol=1.e-9), msg


    def test_function2raster(self):

        if not gdal_installed:
            print 'Warning: GDAL not installed. Not testing '\
                'generic_grid.function2raster'
            return
        
        test_data={'xmin': 0, 'xmax': 5,
                   'ymin': -5, 'ymax':3, 
                   'dx': 1, 'dy': 1,
                   'xx': numpy.array([[0,1,2,3,4],[0,1,2,3,4],[0,1,2,3,4],
                                         [0,1,2,3,4],[0,1,2,3,4],[0,1,2,3,4],
                                         [0,1,2,3,4],[0,1,2,3,4]]),
                   'yy': numpy.array([[-5,-5,-5,-5,-5],[-4,-4,-4,-4,-4],
                                         [-3,-3,-3,-3,-3],[-2,-2,-2,-2,-2],
                                         [-1,-1,-1,-1,-1],[0,0,0,0,0],
                                         [1,1,1,1,1],[2,2,2,2,2]]),
                   'z': numpy.array([[  5. ,   6. ,   9. ,  14. ,  21. ],
                                     [  2.5,   3.5,   6.5,  11.5,  18.5],
                                     [  0. ,   1. ,   4. ,   9. ,  16. ],
                                     [ -2.5,  -1.5,   1.5,   6.5,  13.5],
                                     [ -5. ,  -4. ,  -1. ,   4. ,  11. ],
                                     [ -7.5,  -6.5,  -3.5,   1.5,   8.5],
                                     [-10. ,  -9. ,  -6. ,  -1. ,   6. ],
                                     [-12.5, -11.5,  -8.5,  -3.5,   3.5]])}

        raster_name = 'test.tif'
        generic_grid.function2raster(test_function, filepath=raster_name,
                                     xmin=test_data['xmin'], xmax=test_data['xmax'],
                                     ymin=test_data['ymin'], ymax=test_data['ymax'],
                                     dx=test_data['dx'], dy=test_data['dy'],
                                     format='Gtiff', EPSG_code = 4326, 
                                     proj4string=None, creation_options=[],
                                     gdal_datatype=gdal.GDT_Float32)

        # Now load created raster and get the data to check against test data
        dataset = gdal.Open(raster_name)
        raster_array = numpy.array(dataset.GetRasterBand(1).ReadAsArray())
        msg = 'Generated raster different to expected raster'
        assert numpy.allclose(test_data['z'], raster_array, rtol=1.e-9), msg
        # Remove test file
        os.remove(raster_name)

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_generic_grid,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)




                   
