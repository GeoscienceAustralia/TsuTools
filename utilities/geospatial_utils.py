"""Generic geospatial functionality to be used in TsuTools

Jonathan Griffin, Gareth Davies
Geoscience Australia
March 2015
"""

import numpy

# Check that GDAL interface is installed on import
try:
    import gdal
    import osr
    import ogr
except ImportError, e:
    msg='Failed to import gdal/ogr modules --'\
    + 'perhaps gdal python interface is not installed.'
    raise ImportError, msg

def reproject_point(s_EPSG, t_EPSG, point):
    """
    Function to reproject an individual point into another
    coordinate system. 
    Returns reproject x and y coordinates for the point

    :param s_EPSG:
        EPSG code for source reference system
        See http://spatialreference.org/ref/epsg/
    :param t_EPSG:
        EPSG code for target reference system
        See http://spatialreference.org/ref/epsg/
    :param point:
        List [x,y] of point to be reprojected
    """

    # Create gdal transform object to convert between coordinates systems
    geog = osr.SpatialReference()
    geog.ImportFromEPSG(s_EPSG)
    proj = osr.SpatialReference()
    proj.ImportFromEPSG(t_EPSG)
    transform = osr.CoordinateTransformation(geog, proj)
    point = ogr.CreateGeometryFromWkt("POINT (%.10f %.10f)" % \
                                          (point[0], point[1]))
    point.Transform(transform)
    point_coords = list(point.GetPoint_2D())
    return point_coords[0], point_coords[1]


def make_gdal_grid(data, lons, lats, filepath, format='GTiff', 
              EPSG_code=None, proj4string=None, 
              creation_options=[], gdal_datatype=gdal.GDT_Float32):
    """
    Use GDAL tools to write data to a raster file.
    One of EPSG_code of proj4string must be given to define 
    projection. 

    :param data:
        2D numpy array of z data for raster cell values
    :param lons:
        1D array of longitude or x values
    :param lats:
        1D array of latitude or y values
    :param filepath:
        Path to output raster including file extension
    :param format:
        GDAL raster format code. See 
        http://www.gdal.org/formats_list.html for codes
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

    # Calculate resolution
    xres = lons[1] - lons[0]
    yres = lats[1] - lats[0]
    # Calculate size
    ysize = len(lats)
    xsize = len(lons)

    # Assume data/lons/lats refer to cell centres, 
    # and compute upper left coordinate
    ulx = lons[0] - (xres / 2.)
    uly = lats[lats.shape[0]-1] + (yres / 2.)
    # Need to flip data as writes from top left to bottom right, but
    # data starts from bottom left to top right.
    data = numpy.flipud(data)

    # GDAL magic to make the tif
    driver = gdal.GetDriverByName(format)
    ds = driver.Create(filepath, xsize, ysize, 1, gdal_datatype, 
                       creation_options)

    # Check projection information has been supplied correctly
    srs = osr.SpatialReference()
    if(proj4string is not None):
        srs.ImportFromProj4(proj4string)
    elif(EPSG_code is not None):
        srs.ImportFromEPSG(EPSG_code)
    else:
        raise Exception, 'No spatial reference information given'

    ds.SetProjection(srs.ExportToWkt())

    gt = [ulx, xres, 0, uly, 0, -yres ]
    ds.SetGeoTransform(gt)

    outband = ds.GetRasterBand(1)
    outband.SetNoDataValue(numpy.nan)
    outband.WriteArray(data)

    return
