"""

Authors: Gareth Davies, Geoscience Australia, 2014
         ...
         ...

"""

import os
import shutil
import numpy
import gdal
import struct


def raster_values_at_points(xy, raster_file, band=1, nodata_rel_tol=1.0e-08):
    """Get raster values at point locations.

        INPUT:
        @param xy = numpy array with point locations
        @param raster_file = Filename of the gdal-compatible raster
        @param band = band of the raster to get
        @param nodata_rel_tol = raster pixels are interpreted as nodata if

            abs(raster_values - nodataval) < nodata_rel_tol*abs(nodataval)  

            This approach can compensate for minor rounding errors in
            stored nodata values or nodata pixel values, which I've seen
            a few times. But it could also cause 'real' values to be
            set to nodata if they were extremely close to the nodata value

        OUTPUT:
        1d numpy array with raster values at xy

        This is based on a similar routine
        in ANUGA, also authored by Gareth Davies 2014
    """

    # Raster info
    raster = gdal.Open(raster_file)
    raster_band = raster.GetRasterBand(band)
    raster_band_type = gdal.GetDataTypeName(raster_band.DataType)

    # Projection info
    transform = raster.GetGeoTransform()
    x_origin = transform[0]
    y_origin = transform[3]
    pixel_width = transform[1]
    pixel_height = transform[5]  # Negative

    # Get coordinates in pixel values
    px = ((xy[:, 0] - x_origin) / pixel_width).astype(int)  # x
    py = ((xy[:, 1] - y_origin) / pixel_height).astype(int)  # y

    # Store pixel values
    raster_values = px * 0.

    # Get the right character for struct.unpack
    if (raster_band_type == 'Int16'):
        ctype_name = 'h'
    elif (raster_band_type == 'Float32'):
        ctype_name = 'f'
    elif (raster_band_type == 'Byte'):
        ctype_name = 'B'
    elif (raster_band_type == 'Int32'):
        ctype_name = 'i'
    else:
        print 'unrecognized DataType:', gdal.GetDataTypeName(band.DataType)
        print 'You might need to edit this code to read the data type'
        raise Exception('Stopping')

    # Upper bounds for pixel values, so we can fail gracefully
    xMax = raster.RasterXSize
    yMax = raster.RasterYSize
    if(px.max() < xMax and px.min() >= 0 and
       py.max() < yMax and py.min() >= 0):
        pass
    else:
        msg = 'Trying to extract point values that exceed the raster extent'
        raise Exception(msg)

    # Get values -- seems we have to loop, but it is efficient enough
    for i in range(len(px)):
        xc = int(px[i])
        yc = int(py[i])
        structval = raster_band.ReadRaster(
            xc, yc, 1, 1, buf_type=raster_band.DataType)
        raster_values[i] = struct.unpack(ctype_name, structval)[0]

    # Deal with nodata
    nodataval = raster_band.GetNoDataValue()
    rel_tol = ( abs(raster_values - nodataval) < nodata_rel_tol*abs(nodataval) ) 
    missing = rel_tol.nonzero()[0]
    if len(missing) > 0:
        raster_values[missing] = numpy.nan

    return raster_values

###############################################################################


def make_bathtub_maps(stage_raster_file,
                      elevation_raster_file,
                      ocean_polygon_file,
                      ocean_polygon_file_layername,
                      output_dir='./bathtub_outputs',
                      quiet=False,
                      fill_max_iterations=2000,
                      keep_temp_files=False):
    """Make a 'bathtub' flood map using input stage & elevation rasters +
    an 'ocean polygon'

    It assumes the gdal commandline tools are available via system calls

    @param stage_raster_file
        Filename of raster with stage values
    @param elevation_raster_file
        Filename of raster with elevation values
    @param ocean_polygon_file
        Polygon shapefile. Stage is clipped to this before dilating.
        For tsunami applications it would typically cover ocean areas
    @param ocean_polygon_file_layername
        Layer name (according to ogr) for polygon file. For a shapefile this
        will be the file basename without .shp
    @param output_dir
        Location for output files. Created if necessary.
    @param fill_max_iterations
        The clipped data is flood-filled by at most this many pixels.
        Increase for large rasters
    @param keep_temp_files
        Don't clean up temporary files
    @param quiet
        Don't print commands passed to os.system

    """

    # Make the output directory cleanly
    try:
        os.mkdir(output_dir)
    except:
        pass

    # Preliminary filename definition / creation
    raw_stage_copy = os.path.join(output_dir, 'raw_stage_copy.tif')
    raw_elevation_copy = os.path.join(output_dir, 'raw_elevation_copy.tif')
    poly_mask = os.path.join(output_dir, 'poly_mask.tif')
    stage_clipped = os.path.join(output_dir, 'stage_clipped.tif')
    wet_area = os.path.join(output_dir, 'wet_area.tif')
    stage_flood_filled = os.path.join(
        output_dir, 'stage_flood_filled.tif')
    depth_flood_filled = os.path.join(
        output_dir, 'depth_flood_filled.tif')

    temp_files = [raw_stage_copy, raw_elevation_copy,
                  poly_mask, stage_clipped, wet_area]

    shutil.copyfile(stage_raster_file, raw_stage_copy)
    shutil.copyfile(elevation_raster_file, raw_elevation_copy)

    # Make a 'zero' raster with the same extent as the input rasters
    # Use various tricks to make NA values also 0.0
    zero_raster_cmd = 'gdal_calc.py ' + '-A ' + elevation_raster_file +\
                      ' --calc="((A==A) + (A!=A))*0.0"' +\
                      ' --outfile=' + poly_mask
    if not quiet:
        print zero_raster_cmd
    os.system(zero_raster_cmd)

    # Make a 1/0 mask defining the ocean_polygon as a raster
    polygon_burnin_cmd = 'gdal_rasterize ' + ' -at ' + ' -l ' +\
                         ocean_polygon_file_layername +\
                         ' -burn 1.0 ' +\
                         ocean_polygon_file + ' ' + poly_mask
    if not quiet:
        print polygon_burnin_cmd
    os.system(polygon_burnin_cmd)

    shutil.copyfile(raw_stage_copy, stage_clipped)

    # Dilate the clipped_stage file
    dilate_stage_cmd = 'gdal_fillnodata.py ' +\
        '-md ' + str(int(fill_max_iterations)) + ' ' +\
        stage_clipped + ' -mask ' + poly_mask
    if not quiet:
        print dilate_stage_cmd
    os.system(dilate_stage_cmd)

    # Compute the 'wet-area' for stage/depth calculations
    wet_area_cmd = 'gdal_calc.py ' + ' -A ' + stage_clipped +\
                   ' -B ' + elevation_raster_file + ' --calc="(B<A)"' +\
                   ' --outfile ' + wet_area
    if not quiet:
        print wet_area_cmd
    os.system(wet_area_cmd)

    # Set the filled stage
    filled_stage_cmd = 'gdal_calc.py ' + ' -A ' + stage_clipped +\
        ' -B ' + wet_area +\
        ' -C ' + elevation_raster_file +\
        ' --calc="A*(B==1) + C*(B!=1) "' +\
        ' --outfile ' + stage_flood_filled
    if not quiet:
        print filled_stage_cmd
    os.system(filled_stage_cmd)

    # Convert to depth
    filled_depth_cmd = 'gdal_calc.py ' + ' -A ' + stage_flood_filled +\
        ' -B ' + elevation_raster_file +\
        ' -C ' + wet_area +\
        ' --calc="(A-B)*(C==1) + 0.*( (C==0)+(C!=C))"' +\
        ' --outfile ' + depth_flood_filled
    if not quiet:
        print filled_depth_cmd
    os.system(filled_depth_cmd)

    if not keep_temp_files:
        for filename in temp_files:
            os.remove(filename)

    return


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Make bathtub inundation rasters by clipping stage to a' +
                    ' polygon area, and then dilating the clipped raster.' +
                    ' Creates bathtub stage and bathtub depth rasters')
    parser.add_argument('-stage_raster_file', type=str, default=None,
                        help='Filename of raster with stage values')
    parser.add_argument('-elevation_raster_file', type=str, default=None,
                        help='Filename of raster with elevation values')
    parser.add_argument('-ocean_polygon_file', type=str, default=None,
                        help='Polygon shapefile. Stage is clipped to this ' +
                        ' before dilating. For tsunami applications it' +
                        ' would typically cover ocean areas')
    parser.add_argument(
        '-output_dir', type=str, default='./bathtub_outputs',
        help='Location for output files. Created if necessary.')
    parser.add_argument('-fill_max_iterations', type=int, default=2000,
                        help='The clipped data is flood-filled by at most ' +
                             'this many pixels. Increase for large rasters')
    parser.add_argument('-keep_temp_files', action='store_true', default=False,
                        help=" Don't clean up temporary files")
    parser.add_argument('-quiet', action='store_true', default=False,
                        help="Don't Print commands passed to os.system")

    args = parser.parse_args()

    try:
        ocean_polygon_file_layername = os.path.splitext(os.path.basename(
            args.ocean_polygon_file))[0]

        make_bathtub_maps(args.stage_raster_file,
                          args.elevation_raster_file,
                          args.ocean_polygon_file,
                          ocean_polygon_file_layername,
                          args.output_dir,
                          args.quiet,
                          args.fill_max_iterations,
                          args.keep_temp_files)
    except:
        parser.print_help()
