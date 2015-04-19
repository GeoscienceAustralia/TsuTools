#!/usr/bin/env python

import unittest
import numpy
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as pyplot

import TsuTools
import TsuTools.utilities.generate_event_deformation as generate_event_deformation


class Test_generate_event_deformation(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_distance_points_sphere(self):

        point = numpy.array([104.8, -7.8])

        points = numpy.array([[113.376, -10.625],
                              [113.555, -10.215],
                              [113.725, -9.813],
                              [112.497, -10.410],
                              [112.690, -10.007],
                              [112.872, -9.610],
                              [111.598, -10.281],
                              [111.798, -9.880],
                              [111.992, -9.484],
                              [110.699, -10.205],
                              [110.918, -9.815],
                              [111.132, -9.434],
                              [109.799, -10.105],
                              [110.031, -9.723],
                              [110.262, -9.354],
                              [108.898, -9.999],
                              [109.139, -9.621],
                              [109.377, -9.249],
                              [108.012, -9.842],
                              [108.269, -9.475],
                              [108.515, -9.103],
                              [107.158, -9.539],
                              [107.415, -9.171],
                              [107.667, -8.804],
                              [107.906, -8.434],
                              [106.361, -9.123],
                              [106.630, -8.765],
                              [106.893, -8.405],
                              [105.584, -8.664],
                              [105.859, -8.310],
                              [106.130, -7.957],
                              [104.828, -8.175],
                              [105.109, -7.823],
                              [105.387, -7.476]])

        # For test we computed the answers using R, with the package geosphere,
        # function distHaversine
        distances_computed_using_R = numpy.array(
            [
                992.21893,
                998.20537,
                1005.85919,
                893.43881,
                900.75708,
                909.71254,
                795.78795,
                802.87688,
                812.46502,
                700.82880,
                708.58113,
                719.42673,
                605.92516,
                613.32130,
                624.89564,
                512.28122,
                518.09335,
                528.46046,
                419.64696,
                424.39749,
                433.52039,
                323.37445,
                325.49182,
                334.62249,
                349.10100,
                226.08898,
                228.16760,
                240.02653,
                129.12688,
                129.65313,
                147.52976,
                41.81193,
                34.13634,
                74.04754])

        distances = generate_event_deformation.distance_points_sphere(points,
                                                                      point)

        assert all(abs(distances - distances_computed_using_R) < 1.0e-05)

        # Now test a case where point is in points
        point = numpy.array([105.109, -7.823])

        distances = generate_event_deformation.distance_points_sphere(points,
                                                                      point)

        assert distances.min() == 0.0

        # Check another case with lon > 180 and positive lat
        point[0] = point[0] + 90.0
        point[1] = -point[1]
        points[:, 0] = points[:, 0] + 90.0
        points[:, 1] = -points[:, 1]

        distances_new = generate_event_deformation.distance_points_sphere(
            points,
            point)
        assert all(abs(distances - distances_new) < 1.0e-05)

        # Check another cas with lon < 0
        point[0] = point[0] - 270.0
        points[:, 0] = points[:, 0] - 270.0

        distances_new = generate_event_deformation.distance_points_sphere(
            points,
            point)
        assert all(abs(distances - distances_new) < 1.0e-05)

        return

    def test_nearest_point(self):

        point = numpy.array([105.109, -7.823])

        points = numpy.array([[113.376, -10.625],
                              [113.555, -10.215],
                              [113.725, -9.813],
                              [112.497, -10.410],
                              [112.690, -10.007],
                              [112.872, -9.610],
                              [111.598, -10.281],
                              [111.798, -9.880],
                              [111.992, -9.484],
                              [110.699, -10.205],
                              [110.918, -9.815],
                              [111.132, -9.434],
                              [109.799, -10.105],
                              [110.031, -9.723],
                              [110.262, -9.354],
                              [108.898, -9.999],
                              [109.139, -9.621],
                              [109.377, -9.249],
                              [108.012, -9.842],
                              [108.269, -9.475],
                              [108.515, -9.103],
                              [107.158, -9.539],
                              [107.415, -9.171],
                              [107.667, -8.804],
                              [107.906, -8.434],
                              [106.361, -9.123],
                              [106.630, -8.765],
                              [106.893, -8.405],
                              [105.584, -8.664],
                              [105.859, -8.310],
                              [106.130, -7.957],
                              [104.828, -8.175],
                              [105.109, -7.823],
                              [105.387, -7.476]])

        p1, index = generate_event_deformation.nearest_point_sphere(
            point, points.transpose())

        assert all(abs(point - p1) == 0.0)
        assert index == len(points[:, 0]) - 2

    def test_nearest_point_sphere(self):

        point = numpy.array([105.109, -7.823])

        points = numpy.array([[113.376, -10.625],
                              [113.555, -10.215],
                              [113.725, -9.813],
                              [112.497, -10.410],
                              [112.690, -10.007],
                              [112.872, -9.610],
                              [111.598, -10.281],
                              [111.798, -9.880],
                              [111.992, -9.484],
                              [110.699, -10.205],
                              [110.918, -9.815],
                              [111.132, -9.434],
                              [109.799, -10.105],
                              [110.031, -9.723],
                              [110.262, -9.354],
                              [108.898, -9.999],
                              [109.139, -9.621],
                              [109.377, -9.249],
                              [108.012, -9.842],
                              [108.269, -9.475],
                              [108.515, -9.103],
                              [107.158, -9.539],
                              [107.415, -9.171],
                              [107.667, -8.804],
                              [107.906, -8.434],
                              [106.361, -9.123],
                              [106.630, -8.765],
                              [106.893, -8.405],
                              [105.584, -8.664],
                              [105.859, -8.310],
                              [106.130, -7.957],
                              [104.828, -8.175],
                              [105.109, -7.823],
                              [105.387, -7.476]])

        p1, index = generate_event_deformation.nearest_point_sphere(
            point, points.transpose())

        assert all(abs(point - p1) == 0.0)
        assert index == len(points[:, 0]) - 2

    def test_get_subfaults_in_envelope(self):

        centroid_strike_number = 3
        centroid_dip_number = 2

        strike_numbers = numpy.array([1,
                                      1,
                                      1,
                                      2,
                                      2,
                                      2,
                                      3,
                                      3,
                                      3,
                                      4,
                                      4,
                                      4,
                                      5,
                                      5,
                                      5,
                                      6,
                                      6,
                                      6,
                                      7,
                                      7,
                                      7,
                                      8,
                                      8,
                                      8,
                                      8,
                                      9,
                                      9,
                                      9,
                                      10,
                                      10,
                                      10,
                                      11,
                                      11,
                                      11])

        dip_numbers = numpy.array([1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   4,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3])

        strike_envelope = numpy.array([-5, 5])
        dip_envelope = numpy.array([-1, 1])

        subfaults_inside = generate_event_deformation.get_subfaults_in_envelope(
            centroid_strike_number,
            centroid_dip_number,
            strike_numbers,
            dip_numbers,
            strike_envelope,
            dip_envelope)

        dstrike = strike_numbers[subfaults_inside] - centroid_strike_number
        ddip = dip_numbers[subfaults_inside] - centroid_dip_number
        assert ((dstrike >= strike_envelope[0]) *
                (dstrike <= strike_envelope[1])).all()
        assert ((ddip >= dip_envelope[0]) *
                (ddip <= dip_envelope[1])).all()

    def test_find_subfaults_to_include(self):

        centroid_index = 8

        strike_numbers = numpy.array([1,
                                      1,
                                      1,
                                      2,
                                      2,
                                      2,
                                      3,
                                      3,
                                      3,
                                      4,
                                      4,
                                      4,
                                      5,
                                      5,
                                      5,
                                      6,
                                      6,
                                      6,
                                      7,
                                      7,
                                      7,
                                      8,
                                      8,
                                      8,
                                      8,
                                      9,
                                      9,
                                      9,
                                      10,
                                      10,
                                      10,
                                      11,
                                      11,
                                      11])

        dip_numbers = numpy.array([1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   4,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3])

        centroid_strike_number = strike_numbers[centroid_index]
        centroid_dip_number = dip_numbers[centroid_index]

        # Case 1 -- include everything
        num_subfaults_length = 11
        num_subfaults_width = 4

        subfaults = generate_event_deformation.find_subfaults_to_include(
            strike_numbers, dip_numbers, centroid_index, num_subfaults_length,
            num_subfaults_width)

        assert all(strike_numbers[subfaults] == strike_numbers)
        assert all(dip_numbers[subfaults] == dip_numbers)

        # Case 2 -- include just the centroid
        num_subfaults_length = 1
        num_subfaults_width = 1

        subfaults = generate_event_deformation.find_subfaults_to_include(
            strike_numbers, dip_numbers, centroid_index, num_subfaults_length,
            num_subfaults_width)

        assert strike_numbers[subfaults] == centroid_strike_number
        assert dip_numbers[subfaults] == centroid_dip_number

        # Case 3 -- a case where we can't match all the constraints

        num_subfaults_length = 7
        num_subfaults_width = 2

        # Centroid neighbours include missing values
        centroid_index = (dip_numbers == 4).nonzero()[0]

        subfaults = generate_event_deformation.find_subfaults_to_include(
            strike_numbers, dip_numbers, centroid_index, num_subfaults_length,
            num_subfaults_width)

        assert (centroid_index in subfaults.tolist())

        empirical_length = strike_numbers[subfaults].max() -\
            strike_numbers[subfaults].min() + 1
        tol = 0.7 * num_subfaults_length
        assert (abs(empirical_length - num_subfaults_length) <= tol)

        empirical_width = dip_numbers[subfaults].max() -\
            dip_numbers[subfaults].min() + 1
        tol = 0.5 * num_subfaults_width
        assert (abs(empirical_width - num_subfaults_width) <= tol)

        empirical_area = len(subfaults)
        desired_area = num_subfaults_width * num_subfaults_length
        tol = desired_area * 0.2
        assert (abs(empirical_area - desired_area) <= tol)

    def test_plot_subfaults(self):

        strike_numbers = numpy.array([1,
                                      1,
                                      1,
                                      2,
                                      2,
                                      2,
                                      3,
                                      3,
                                      3,
                                      4,
                                      4,
                                      4,
                                      5,
                                      5,
                                      5,
                                      6,
                                      6,
                                      6,
                                      7,
                                      7,
                                      7,
                                      8,
                                      8,
                                      8,
                                      8,
                                      9,
                                      9,
                                      9,
                                      10,
                                      10,
                                      10,
                                      11,
                                      11,
                                      11])

        dip_numbers = numpy.array([1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   4,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3,
                                   1,
                                   2,
                                   3])

        # Let's try making the centroids everywhere
        pp = PdfPages('multipage.pdf')

        for i in range(len(strike_numbers)):
            centroid_index = i

            num_subfaults_length = 6
            num_subfaults_width = 2

            f, axarr = pyplot.subplots(2, 2)
            for j in [0, 1]:
                for k in [0, 1]:
                    subfaults = generate_event_deformation.find_subfaults_to_include(
                        strike_numbers,
                        dip_numbers,
                        centroid_index,
                        num_subfaults_length,
                        num_subfaults_width,
                        prefer_down_dip_expansion=j,
                        prefer_along_strike_expansion=k)
                    # pyplot.subplot(2,2,2*j+k)
                    axarr[
                        j,
                        k].plot(
                        strike_numbers,
                        dip_numbers,
                        'p',
                        color='black')
                    axarr[
                        j,
                        k].plot(
                        strike_numbers[subfaults],
                        dip_numbers[subfaults],
                        'p',
                        color='orange')
                    axarr[
                        j, k].set_xlim(
                        (strike_numbers.min() - 1, strike_numbers.max() + 1))
                    axarr[
                        j, k].set_ylim(
                        (dip_numbers.min() - 1, dip_numbers.max() + 1))
                    axarr[
                        j,
                        k].set_title(
                        'Up Dip: ' +
                        str(j) +
                        ' / Along Strk:' +
                        str(k))
                    axarr[
                        j,
                        k].plot(
                        strike_numbers[centroid_index],
                        dip_numbers[centroid_index],
                        'p',
                        color='red')
            pp.savefig()

        pp.close()

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_generate_event_deformation, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
