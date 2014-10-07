mycmd = "python bathtub.py " +\
        "-stage_raster_file  'test_data/Res_500_initial_stage.tif' " +\
        "-elevation_raster_file 'test_data/Res_500_elevation_less_100.tif' " +\
        "-ocean_polygon_file 'test_data/ocean_polygon/ocean_polygon.shp' " +\
        ""
        #"-keep_temp_files"

import os
os.system(mycmd)
