import os
import shutil

def make_bathtub_maps(stage_raster_file,
                      elevation_raster_file,
                      ocean_polygon_file,
                      ocean_polygon_file_layername,
                      output_dir='./bathtub_outputs',
                      verbose=True,
                      fill_max_iterations=2000,
                      clean_temp_files=True):
    """

    Make a 'bathtub' flood map using input stage & elevation rasters + 
    an 'ocean polygon'

    The computations proceed as:
    1) The stage is first 'masked' to the 'ocean polygon'.
    2) gdal_fillnodata.py is used to expand (dilate) the stage raster values,
       making the 'bathtub stage'.
    3) The 'bathtub-stage' is then set to zero where 'elevation > dilated stage'.
       This approach will not respect 'topographic barriers', so all elevation
       cells < nearest stage will be filled
    4) The 'bathtub depth' is computed from the bathtub stage + input elevation

    It assumes the standard gdal commandline tools are available via system calls

    """

    # Make the output directory cleanly
    try:
        os.mkdir(output_dir)
    except:
        pass

    ###########################################################################
    #
    # Make a raster mask which is 1.0 inside the ocean_polygon
    #
    ###########################################################################

    # Preliminary filename definition / creation
    raw_stage_copy = os.path.join(output_dir, 'raw_stage_copy.tif')
    poly_mask = os.path.join(output_dir, 'poly_mask.tif')
    stage_clipped = os.path.join(output_dir, 'stage_clipped.tif')
    stage_flood_filled_clipped = os.path.join(
        output_dir, 'stage_flood_filled_clipped.tif')
    depth_flood_filled_clipped = os.path.join(
        output_dir, 'depth_flood_filled_clipped.tif')

    temp_files = [raw_stage_copy, poly_mask, stage_clipped]

    shutil.copyfile(stage_raster_file, raw_stage_copy)
    shutil.copyfile(stage_raster_file, poly_mask)

    # Make a 'zero' raster with the same extent as the stage raster
    zero_raster_cmd = 'gdal_calc.py ' + '-A ' + raw_stage_copy +\
                      ' --calc=A*0.' +\
                      ' --outfile=' + poly_mask + ' --NoDataValue=0.0'
    if verbose:
        print zero_raster_cmd
    os.system(zero_raster_cmd)

    # Make a 1/0 mask defining the ocean_polygon as a raster
    polygon_burnin_cmd = 'gdal_rasterize ' + ' -at ' + ' -l ' +\
                         ocean_polygon_file_layername +\
                         ' -burn 1.0 ' +\
                         ocean_polygon_file + ' ' + poly_mask
    if verbose:
        print polygon_burnin_cmd
    os.system(polygon_burnin_cmd)

    # Make the 'clipped stage' file
    clip_stage_cmd = 'gdal_calc.py ' +\
                     ' -A ' + raw_stage_copy +\
                     ' -B ' + poly_mask + ' --calc=A*B ' +\
                     ' --outfile ' + stage_clipped
    if verbose:
        print clip_stage_cmd
    os.system(clip_stage_cmd)

    # Fill the clipped_stage file
    fill_stage_cmd = 'gdal_fillnodata.py ' +\
                     '-md ' + str(int(fill_max_iterations)) + ' ' +\
                     stage_clipped + ' -mask ' + poly_mask
    if verbose:
        print fill_stage_cmd
    os.system(fill_stage_cmd)

    # Clip the filled stage based on the elevation
    clip_stage_cmd = 'gdal_calc.py ' + ' -A ' + stage_clipped +\
                     ' -B ' + elevation_raster_file + ' --calc="A*(B<A)"' +\
                     ' --outfile ' + stage_flood_filled_clipped
    if verbose:
        print clip_stage_cmd
    os.system(clip_stage_cmd)

    # Convert to depth
    depth_cmd = 'gdal_calc.py ' + ' -A ' + stage_flood_filled_clipped +\
                ' -B ' + elevation_raster_file + ' --calc="(A-B)*(B<A)"' +\
                ' --outfile ' + depth_flood_filled_clipped
    if verbose:
        print depth_cmd
    os.system(depth_cmd)

    if clean_temp_files:
        for filename in temp_files:
            os.remove(filename)

    return


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Make bathtub inundation rasters by clipping stage to a' +\
                    ' polygon area, and then dilating the clipped raster')
    parser.add_argument('-stage_raster_file', type=str, default = None,
                        help='Filename of raster with stage values')
    parser.add_argument('-elevation_raster_file', type=str, default=None,
                        help='Filename of raster with elevation values')
    parser.add_argument('-ocean_polygon_file', type=str, default=None, 
                        help='Polygon shapefile. Stage is clipped to this ' +\
                        ' before dilating. For tsunami applications it' +\
                        ' would typically cover ocean areas')
    parser.add_argument(
        '-output_dir', type=str, default='./bathtub_outputs',
        help='Location for output files. Created if necessary.')
    parser.add_argument('-fill_max_iterations', type=int, default=2000,
                        help='The clipped data is flood-filled by at most ' +\
                             'this many pixels. Increase for large rasters')
    parser.add_argument('-clean_temp_files', type=bool, default=True,
                        help='Clean up temporary files')
    parser.add_argument('-verbose', type=bool, default=True,
                        help='Print commands passed to os.system')

    args = parser.parse_args()

    ocean_polygon_file_layername = os.path.splitext(os.path.basename(
        args.ocean_polygon_file))[0]

    make_bathtub_maps(args.stage_raster_file,
                      args.elevation_raster_file,
                      args.ocean_polygon_file,
                      ocean_polygon_file_layername,
                      args.output_dir,
                      args.verbose,
                      args.fill_max_iterations,
                      args.clean_temp_files)

