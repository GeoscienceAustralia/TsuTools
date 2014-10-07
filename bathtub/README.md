bathtub
============

This is a commandline script to make a basic 'bathtub' inundation map from input stage and elevation raster data, and a polygon shapefile defining the extent of the 'ocean' stages which are dilated to 'fill' the other areas.

For us, the purpose is to compute 'bathtub' inundation maps from tsunami models, based only on the
computed 'ocean' stage. 

Bathtub maps are a useful reference case for comparison with detailed inundation maps, as while the latter can be more accurate, they can also be affected by errors due to e.g. mesh resolution, friction definition, onshore elevation data, etc.

The computations proceed as:
----------------------------
* The stage is first 'masked' to the 'ocean polygon'.
* gdal_fillnodata.py is used to expand (dilate) the stage raster values, making the 'bathtub stage'. See gdal_fillnodata for information on how dilation works.
* The 'bathtub-stage' is then set to zero where 'elevation > dilated stage'. ``` This approach will not respect 'topographic barriers', so all elevation cells < nearest stage will be filled ```
* The 'bathtub depth' is computed from the bathtub stage + input elevation

*It assumes the standard gdal commandline tools are available via system calls*

Getting started:
----------------
Run

    python bathtub.py
to see the command-line help, and also look at the example in *example_usage.py*

    python example_usage.py
