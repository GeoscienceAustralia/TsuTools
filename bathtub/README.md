bathtub
============

This is a commandline script to make a basic 'bathtub' inundation map from input stage and elevation raster data, and a polygon shapefile defining the extent of the 'ocean' stages which are dilated to 'fill' the other areas.

For us, the purpose is to compute 'bathtub' inundation maps from tsunami models, based only on the computed 'ocean' stage. 

The maps set the stage in each cell based on their nearest 'ocean' raster value. ``` This approach will not respect 'topographic barriers' -- all elevation cells < nearest stage will be filled ```

Bathtub maps are a useful reference case for comparison with inundation maps based on hydrodynamic models, as while the latter can be more accurate, they can also be affected by errors due to e.g. mesh resolution, friction definition, onshore elevation data, etc.

The computations proceed as:
----------------------------
* The stage is first 'masked' to the 'ocean polygon'.
* gdal_fillnodata.py is used to expand (dilate) the 'masked stage' raster values, making the 'bathtub stage'. Basically, this means pixels have a stage set based on their 'nearest' stage in the ocean polygon (irrespective of topographic barriers). See gdal_fillnodata for information on how dilation works. 
* The 'bathtub stage' is finally set to the max of the previously computed bathtub stage, and the input elevation
* The 'bathtub depth' is computed from the difference between the bathtub stage and the input elevation, with a min of 0.0

No-data values:  
---------------

* Nodata pixels in the input stage raster may be flood-filled (if their 'bathtub stage' > elevation). 
* Nodata pixels in the input elevation raster should result in nodata values in the output rasters

Dependencies: 
-------------
We assume the standard gdal commandline tools are available via system calls

Getting started:
----------------
Run

    python bathtub.py
to see the command-line help, and also look at the example in *example_usage.py*

    python example_usage.py

Bugs:
----
Please use github pull requests to contribute bug fixes. Alternatively email Gareth Davies, grothered ( at gmail dot com)

