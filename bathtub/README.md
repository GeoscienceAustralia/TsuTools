bathtub
============

This contains a commandline script to make a basic 'bathtub' inundation map from input stage and elevation raster data, and a polygon shapefile defining the extent of the 'ocean' stages which are dilated to 'fill' the other areas.

The bathtub module also contains a function to extract raster cell values at points, and to run the bathtub computations

For us, the purpose is to compute 'bathtub' inundation maps from tsunami models, based only on the computed 'ocean' stage. Some tsunami inundation studies compute only the offshore elevations, and then extrapolate this onshore with a bathtub type methodology.

The algorithm sets the stage in each cell based on its 'nearby ocean' raster values using gdal_fillnodata.py. ``` This approach will not respect 'topographic barriers```

Bathtub maps are a useful reference case for comparison with inundation maps based on hydrodynamic models, as while the latter can be more accurate, they can also be affected by errors due to e.g. mesh resolution, friction definition, onshore elevation data, etc. Comparison with a bathtub map may help you to discover problems such as the latter, and/or to better understand the inundation flow dynamics.

The computations proceed as:
----------------------------
* The stage is first 'masked' to the 'ocean polygon'.
* gdal_fillnodata.py is used to expand (dilate) the 'masked stage' raster values, making the 'bathtub stage'. Basically, this means pixels have a stage set based on their 'nearby' stages in the ocean polygon (irrespective of topographic barriers). The GDALFillNoData documentation describes the dilation algorithm thus:

    *This algorithm will interpolate values for all designated nodata
    pixels..... For each pixel a four direction conic search is done to find values
    to interpolate from (using inverse distance weighting) ..... This algorithm is
    generally suitable for interpolating missing regions of fairly continuously
    varying rasters (such as elevation models for instance). It is also suitable
    for filling small holes and cracks in more irregularly varying images (like
    airphotos). It is generally not so great for interpolating a raster from sparse
    point data - see the algorithms defined in gdal_grid.h for that case.*
* The 'bathtub stage' is finally set to the max of the previously computed bathtub stage, and the input elevation
* The 'bathtub depth' is computed from the difference between the bathtub stage and the input elevation, with a min of 0.0

No-data values:  
---------------

* Nodata pixels in the input stage raster may be flood-filled (if their 'bathtub stage' > elevation). This is done because some stage rasters use 'nodata' to denote dry areas [the other reasonable approach is to set stage = elevation in dry areas -- this is the approach taken for output tifs in the bathtub routine]
* Nodata pixels in the input elevation raster should result in nodata values in the output rasters

Dependencies: 
-------------
We assume the standard gdal commandline tools are available via system calls, and that numpy, gdal, and standard python library modules are available.

Getting started:
----------------
Run

    python bathtub.py
to see the command-line help, and also look at the example in *example_usage.py*

    python example_usage.py

Run the unit tests with

    python test_bathtub.py


Bugs:
----
Please use github pull requests to contribute bug fixes. Alternatively email Gareth Davies, gareth.davies.ga.code@gmail.com, although we make no guarentees to maintain this software.

Issues: 
--------
* The gdal_fillnodata algorithm might be prone to minor artefacts / strange interpolation? I have not yet seen this of significant magnitude, but strongly suggest a visual comparison of the filled stage raster with the original stage.
