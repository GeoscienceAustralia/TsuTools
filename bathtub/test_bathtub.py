#!/usr/bin/env python

import unittest
import bathtub
import numpy


class Test_bathtub(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_raster_values_at_points(self):
        """ Test raster_values_at_points works ok

        """
        # These point values were manually extracted (to within the reported
        # decimal places) in GIS
        raster_file = './test_data/Res_500_initial_stage.tif'
        point_values = [[445496, 9387281, 0.96923],
                        [352454, 9389770, -0.02730],
                        [344754, 9401707, numpy.nan]]
        point_values = numpy.array(point_values)

        extracted_values = bathtub.raster_values_at_points(
            point_values[:, 0:2], raster_file)

        assert abs(extracted_values[0] - point_values[0, 2]) < 1.0e-04
        assert abs(extracted_values[1] - point_values[1, 2]) < 1.0e-04
        # To check for nan, note that nan != nan
        assert extracted_values[2] != extracted_values[2]

        return

    def test_example_usage(self):
        """Check we can run example_usage without errors,
           and that some extracted values behave as expected

        """
        cmd = 'python example_usage.py'
        import os
        os.system(cmd)

        stage_filled_file = './bathtub_outputs/stage_flood_filled.tif'
        depth_filled_file = './bathtub_outputs/depth_flood_filled.tif'
        elevation_file = './test_data/Res_500_elevation_less_100.tif'
        stage_file = './test_data/Res_500_initial_stage.tif'

        # Make points inside the ocean
        inside_ocean_points = [[462900., 9325505.],
                               [491114., 9458604.],
                               [654703., 9349545.]]
        inside_ocean_points = numpy.array(inside_ocean_points)

        # Inside the ocean, stage_filled should = original stage
        filled_stage_ocean = bathtub.raster_values_at_points(
            inside_ocean_points, stage_filled_file)
        raw_stage_ocean = bathtub.raster_values_at_points(
            inside_ocean_points, stage_file)
        assert numpy.allclose(filled_stage_ocean, raw_stage_ocean)

        # Check that ocean has elevation+depth=stage
        filled_depth_ocean = bathtub.raster_values_at_points(
            inside_ocean_points, depth_filled_file)
        elevation_ocean = bathtub.raster_values_at_points(
            inside_ocean_points, elevation_file)
        # Minor (rounding) error does appear, probably due to
        # the residual being a small fraction of the depth/elevation
        # from which it was derived
        assert numpy.allclose(filled_depth_ocean + elevation_ocean,
                              raw_stage_ocean, atol=1.0e-03)

        # Check that the bathtub stage values are 'fairly similar'
        # to the nearby ocean values. I am not working from a rigorous
        # understanding of gdal_fillnodata's algorithm [although the
        # source is available]. Anyway we can check for general similarity
        nearby_pairs_inside = numpy.array(
            [[670742., 9344373.], [392230., 9381899.]])
        nearby_pairs_outside = numpy.array(
            [[687160., 9344373.], [389299., 9388348.]])

        n_in = bathtub.raster_values_at_points(
            nearby_pairs_inside, stage_file)
        n_out = bathtub.raster_values_at_points(
            nearby_pairs_inside, stage_filled_file)

        assert numpy.allclose(n_in, n_out, atol=1.0e-02)
        return


if __name__ == "__main__":
    suite = unittest.makeSuite(Test_bathtub, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
