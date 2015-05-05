#!/usr/bin/env python

import unittest
from strasser2010 import *

class Test_strasser2010(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_calculate_dimensions(self):
        """ Test that correct dimensions are returned for range of magnitudes
        """

        test_data = {'magnitude': [7.0, 8.0, 9.1],
                     'length': [41.4954042634, 159.5879147237, 702.2633437358],
                     'width': [37.5837404288, 84.3334757764, 205.1634530751],
                     'area': [1541.7004529496, 13803.8426460288, 
                            153886.3149819950],
                     'sigma_length_plusone': [62.80583588, 241.5460834, \
                                                  1062.918583],
                     'sigma_width_plusone': [55.97576015, 125.6029964, \
                                                 305.5624616],
                     'sigma_area_plusone': [3104.559588, 27797.13268, \
                                                309884.6042],
                     'sigma_length_minusone': [27.41574172, 105.4386896, \
                                                   463.98079],
                     'sigma_width_minusone': [25.23480772, 56.6239289, \
                                                  137.752662],
                     'sigma_area_minusone': [765.5966069, 6854.882265, \
                                                 76418.7624]}

        i = 0
        for magnitude in test_data['magnitude']:
            length, width, area, sigma_length_plusone, sigma_width_plusone, \
            sigma_area_plusone, sigma_length_minusone, sigma_width_minusone, \
            sigma_area_minusone = calculate_dimensions(magnitude)
            msg = 'Expected value of %.10f, got %.10f for length for \
                   magnitude %.1f' % (test_data['length'][i], length, magnitude)
            assert numpy.allclose(test_data['length'][i], length, rtol=1.e-9),\
                   msg
            msg = 'Expected value of %.10f, got %.10f for width for \
                   magnitude %.1f' % (test_data['width'][i], width, magnitude)
            assert numpy.allclose(test_data['width'][i], width, rtol=1.e-9),\
                   msg
            msg = 'Expected value of %.10f, got %.10f for area for \
                   magnitude %.1f' % (test_data['area'][i], area, magnitude)
            assert numpy.allclose(test_data['area'][i], area, rtol=1.e-9), msg

            msg = 'Expected value of %.10f, got %.10f for length + 1 sigma for \
                   magnitude %.1f' % (test_data['sigma_length_plusone'][i], \
                                          sigma_length_plusone, magnitude)
            assert numpy.allclose(test_data['sigma_length_plusone'][i], \
                                      sigma_length_plusone, rtol=1.e-9), msg
            msg = 'Expected value of %.10f, got %.10f for width + 1 sigma for \
                   magnitude %.1f' % (test_data['sigma_width_plusone'][i], \
                                          sigma_width_plusone, magnitude)
            assert numpy.allclose(test_data['sigma_width_plusone'][i], \
                                      sigma_width_plusone, rtol=1.e-9), msg
            msg = 'Expected value of %.10f, got %.10f for area + 1 sigma for \
                   magnitude %.1f' % (test_data['sigma_area_plusone'][i], \
                                          sigma_area_plusone, magnitude)
            assert numpy.allclose(test_data['sigma_area_plusone'][i], \
                                      sigma_area_plusone, rtol=1.e-9), msg

            msg = 'Expected value of %.10f, got %.10f for length - 1 sigma for \
                   magnitude %.1f' % (test_data['sigma_length_minusone'][i], \
                                          sigma_length_minusone, magnitude)
            assert numpy.allclose(test_data['sigma_length_minusone'][i], \
                                      sigma_length_minusone, rtol=1.e-9), msg
            msg = 'Expected value of %.10f, got %.10f for width - 1 sigma for \
                   magnitude %.1f' % (test_data['sigma_width_plusone'][i], \
                                          sigma_width_minusone, magnitude)
            assert numpy.allclose(test_data['sigma_width_minusone'][i], \
                                      sigma_width_minusone, rtol=1.e-9), msg
            msg = 'Expected value of %.10f, got %.10f for area - 1 sigma for \
                   magnitude %.1f' % (test_data['sigma_area_minusone'][i], \
                                          sigma_area_minusone, magnitude)
            assert numpy.allclose(test_data['sigma_area_minusone'][i], \
                                      sigma_area_minusone, rtol=1.e-9), msg

            i+=1

    def test_calculate_magnitude(self):
        """ Test that correct magnitudes are returned for given dimensions
        """

        test_length_data = {'length': [50.0, 160., 700.0],
                            'magnitude': [7.2329662460, 7.9361350159, 
                                          8.8283764717],
                            'magnitude_plusone_sigma': [7.509966246,
                                                        8.213135016,
                                                        9.105376472],
                            'magnitude_minusone_sigma': [6.955966246,
                                                         7.659135016,
                                                         8.551376472]}
        test_width_data = {'width': [40.0, 80.0, 200.0],
                           'magnitude': [7.3017182843, 7.8450774265, 
                                         8.5633591422],
                           'magnitude_plusone_sigma': [7.693718284,
                                                       8.237077427,
                                                       8.955359142],
                           'magnitude_minusone_sigma': [6.909718284,
                                                        7.453077427,
                                                        8.171359142]}
        test_area_data = {'area': [1500.0, 14000.0, 150000.0],
                          'magnitude': [7.1279732052, 7.9486243182, 
                                        8.8199732052],
                          'magnitude_plusone_sigma': [7.413973205,
                                                      8.234624318,
                                                      9.105973205],
                          'magnitude_minusone_sigma': [6.841973205,
                                                       7.662624318,
                                                       8.533973205]}
        i = 0
        for length in test_length_data['length']:
            magnitude,magnitude_plusone_sigma, magnitude_minusone_sigma \
                = calculate_magnitude('length', length)
            msg = 'Expected magnitude value of %.10f, got %.10f for \
                   length %.1f' % (test_length_data['magnitude'][i], magnitude,
                                   length)
            assert numpy.allclose(test_length_data['magnitude'][i], 
                                  magnitude, rtol=1.e-9), msg
            msg = 'Expected magnitude + one sigma value of %.10f, got %.10f' \
                'for length %.1f' % \
                (test_length_data['magnitude_plusone_sigma'][i], \
                     magnitude_plusone_sigma,length)
            assert numpy.allclose(test_length_data['magnitude_plusone_sigma'][i], 
                                  magnitude_plusone_sigma, rtol=1.e-9), msg
            msg = 'Expected magnitude - one sigma value of %.10f, got %.10f' \
                'for length %.1f' % \
                (test_length_data['magnitude_minusone_sigma'][i], \
                     magnitude_minusone_sigma,length)
            assert numpy.allclose(test_length_data['magnitude_minusone_sigma'][i], 
                                  magnitude_minusone_sigma, rtol=1.e-9), msg
            i+=1

        i = 0
        for width in test_width_data['width']:
            magnitude, magnitude_plusone_sigma, magnitude_minusone_sigma = \
                calculate_magnitude('width', width)
            msg = 'Expected magnitude value of %.10f, got %.10f for \
                   width %.1f' % (test_width_data['magnitude'][i], magnitude,
                                   width)
            assert numpy.allclose(test_width_data['magnitude'][i], 
                                  magnitude, rtol=1.e-9), msg
            msg = 'Expected magnitude + one sigma value of %.10f, got %.10f ' \
                'for width %.1f' % (test_width_data['magnitude_plusone_sigma'][i],\
                                        magnitude_plusone_sigma,width)
            assert numpy.allclose(test_width_data['magnitude_plusone_sigma'][i], 
                                  magnitude_plusone_sigma, rtol=1.e-9), msg
            msg = 'Expected magnitude - one sigma value of %.10f, got %.10f ' \
                'for width %.1f' % (test_width_data['magnitude_minusone_sigma'][i],\
                                        magnitude_minusone_sigma,width)
            assert numpy.allclose(test_width_data['magnitude_minusone_sigma'][i], 
                                  magnitude_minusone_sigma, rtol=1.e-9), msg
            i+=1

        i = 0
        for area in test_area_data['area']:
            magnitude, magnitude_plusone_sigma, magnitude_minusone_sigma  = \
                calculate_magnitude('area', area)
            msg = 'Expected magnitude value of %.10f, got %.10f for \
                   area %.1f' % (test_area_data['magnitude'][i], magnitude,
                                   area)
            assert numpy.allclose(test_area_data['magnitude'][i], 
                                  magnitude, rtol=1.e-9), msg
            msg = 'Expected magnitude + one sigma value of %.10f, got %.10f ' \
                'for area %.1f' % (test_area_data['magnitude_plusone_sigma'][i], \
                                       magnitude_plusone_sigma, area)
            assert numpy.allclose(test_area_data['magnitude_plusone_sigma'][i], 
                                  magnitude_plusone_sigma, rtol=1.e-9), msg
            msg = 'Expected magnitude - one sigma value of %.10f, got %.10f ' \
                'for area %.1f' % (test_area_data['magnitude_minusone_sigma'][i], \
                                       magnitude_minusone_sigma, area)
            assert numpy.allclose(test_area_data['magnitude_minusone_sigma'][i], 
                                  magnitude_minusone_sigma, rtol=1.e-9), msg
            i+=1



if __name__ == "__main__":
    suite = unittest.makeSuite(Test_strasser2010, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
