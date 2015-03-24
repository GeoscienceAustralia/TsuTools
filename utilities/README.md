utilities
==================================
Contains a variety of re-usable scripts

========================================
generate_event_deformation.py

Usage:
python generate_event_deformation.py <source_zone> <Mw> <centroid_longitude> <centroid_latitude>
source_zone must be the name of a folder in the current working directory that contains subfault information in i_invall-<source_zone> and precalculated unit source deformation grd files ./<source_zone>/deformation_data/<source_zone>-<subfault_index>-1m.grd


========================================
geospatial_utils.py

Contains various GIS functions, including functions for writing data to raster files using gdal and functions for projecting points to different coordinate systems