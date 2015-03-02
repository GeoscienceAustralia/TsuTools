"""Create a grid of points over which to evaluate a function
for generating an initial surface deformation

Jonathan Griffin
Geoscience Australia
February 2015
"""

import numpy
from TsuTools.utilities import geospatial_utils

# Check that GDAL interface is installed on import
try:
    import gdal
    import osr
except ImportError, e:
    msg='Failed to import gdal/ogr modules --'\
    + 'perhaps gdal python interface is not installed.'
    raise ImportError, msg

def create_grid(xmin, xmax, ymin, ymax, 
                dx, dy):
    """
    Create a grid of points for evaluating a function

    :param xmin:
        Minimum coordinate along the x-axis
    :param xmax:
        Maximum coordinate along the x-axis
    :param ymin:
        Minimum coordinate along the y-axis
    :param yax:
        Maximum coordinate along the y-axis
    :param xmin:
        Minimum coordinate along the x-axis
    :param dx:
        grid spacing along the x-axis
    :param dy:
        grid spacing along the y-axis
#    :param projection:
#        Geographic projection if using projected coordinates
    """

    x = numpy.arange(xmin, xmax, dx)
    y = numpy.arange(ymin, ymax, dy)
    xx, yy = numpy.meshgrid(x,y)
    
    return x, y, xx, yy

def evaluate_function_over_grid(func,
                                gridx=None, gridy=None,
                                xmin=None, xmax=None,
                                ymin=None, ymax=None,
                                dx=None, dy=None):
    """
    Evaluate a given function over a grid. 
    Either uses existing grid or creates new one
    
    :param func:
        Function F(x,y) to be evaluated
    :param gridx:
        Array of x coordinates of grid
    :param gridy:
        Array of y coordinates of grid
    :param xmin:
        Minimum coordinate along the x-axis
    :param xmax:
        Maximum coordinate along the x-axis
    :param ymin:
        Minimum coordinate along the y-axis
    :param yax:
        Maximum coordinate along the y-axis
    :param xmin:
        Minimum coordinate along the x-axis
    :param dx:
        grid spacing along the x-axis
    :param dy:
        grid spacing along the y-axis

    """

    if gridx is None and gridy is None:
        x, y, gridx, gridy = create_grid(xmin, xmax, ymin, ymax, dx, dy)
    else:
        x = gridx[0,:]
        y = gridy[:,0]
    
    z = func(gridx, gridy)
    return x, y, gridx, gridy, z

def function2raster(func, filepath=None, gridx=None, gridy=None,
                    xmin=None, xmax=None, ymin=None, ymax=None,
                    dx=None, dy=None, format='GTiff',
                    EPSG_code=None, proj4string=None, 
                    creation_options=[], 
                    gdal_datatype=gdal.GDT_Float32):
    """
    Takes a function, evaluates it over an x,y range and creates
    a raster of the result using the GDAL tools. Can either take
    x, y coordinates for a grid (gridx, gridy) or just the range of the axes
    (xmin, xmax, ymin, ymax, dx, dy).
    
    :param func:
        Function F(x,y) to be evaluated
    :param filepath:
        Path to output raster with file extension
    :param gridx:
        Array of x coordinates of grid
    :param gridy:
        Array of y coordinates of grid
    :param xmin:
        Minimum coordinate along the x-axis
    :param xmax:
        Maximum coordinate along the x-axis
    :param ymin:
        Minimum coordinate along the y-axis
    :param yax:
        Maximum coordinate along the y-axis
    :param xmin:
        Minimum coordinate along the x-axis
    :param dx:
        grid spacing along the x-axis
    :param dy:
        grid spacing along the y-axis
    :param EPSG_code:
        Integer code with projection information in EPSG format
        See http://spatialreference.org/ref/epsg/
    :param proj4string:
        proj4string with projection information
        See https://trac.osgeo.org/proj/wiki/GenParms
    :param creation_options:
        List of additional GDAL raster creation options
    :param gdal_datatype:
        Format for output data 
        See http://www.gdal.org/gdal_8h.html#a22e22ce0a55036a96f652765793fb7a4
    """

    if filepath is None:
        filepath='raster.tif'
    
    if gridx is None and gridy is None:
        x, y, gridx, gridy, z = evaluate_function_over_grid(func, 
                                                            xmin=xmin, 
                                                            xmax=xmax,
                                                            ymin=ymin,
                                                            ymax=ymax,
                                                            dx=dx, dy=dy)
    else:
        x, y, gridx, gridy, z = evaluate_function_over_grid(func, 
                                                            gridx=gridx,
                                                            gridy=gridy)
                                    
    geospatial_utils.make_gdal_grid(z, x, y, filepath, format=format,
                   EPSG_code = EPSG_code, proj4string=proj4string,
                   creation_options=[], gdal_datatype=gdal.GDT_Float32)

    return
