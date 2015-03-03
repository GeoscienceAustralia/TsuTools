"""Generate initial sea surface elevation condition for a submarine mass
slide

Jonathan Griffin
Geoscience Australia
March 2015
"""

from subprocess import call

from TsuTools.tsusrc.TsuSource.mass_movement import smf
from TsuTools.tsusrc.TsuSource.grid import generic_grid
# Check that GDAL interface is installed on import
try:
    import gdal
    import osr
except ImportError, e:
    msg='Failed to import gdal/ogr modules --'\
    + 'perhaps gdal python interface is not installed.'
    raise ImportError, msg
from TsuTools.utilities import geospatial_utils

deformation_raster = 'initial_condition.tif'
# Slide parameters, defined by user
length = 600.0
dep = 150.0
theta = 9.0 # slope (degrees)
thk = 15.0
width = 340.0
kappa = 3.0
kappad = 0.8
zsmall = 0.01
alpha=-90.0 # 0.0 means slide is from west to east !!!
gravity=9.8
gamma=1.85
massco=1
dragco=1
frictionco=0
psi=0
x0 = 73.65 # Slide origin (x coordinate)
y0 = -53.04 # Slide origin (y coordinate)

# Parameters for initial conditions raster
xmin=73.0
xmax=74.5
ymin=-53.5
ymax=-52.5
dx=100 # In metres
dy=100 # In metres
geographic = True # If true, will project to UTM to evaluate 
# slide function before reprojecting to geographic coordinates.
# Assumes geographic coordinates are in WGS84 datum
EPSG_code_wgs84 = 4326
UTM_zone = -43
EPSG_code = 32743

###############################################################
# No edits should be required below here

# Transform slide and grid locations to utm
x0, y0 = geospatial_utils.reproject_point(EPSG_code_wgs84, EPSG_code, [x0, y0])

xmin, ymin = geospatial_utils.reproject_point(EPSG_code_wgs84, EPSG_code, \
                                                  [xmin, ymin])
xmax, ymax = geospatial_utils.reproject_point(EPSG_code_wgs84, EPSG_code, \
                                                  [xmax, ymax])

# Generate function
slide_function = smf.slide_tsunami(length=length, depth=dep, slope=theta,
                                   width = width, thickness=thk,
                                   x0=x0, y0=y0, alpha=alpha, 
                                   gravity=gravity, gamma=gamma,
                                   massco=massco, dragco=dragco, 
                                   frictionco=frictionco, psi=psi,
                                   kappa=kappa, kappad=kappad, 
                                   zsmall=zsmall, verbose=False)

# Evaluate function over region of interest
generic_grid.function2raster(slide_function, filepath = deformation_raster,
                             xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                             dx=dx, dy=dy, EPSG_code=EPSG_code)

# Use GDAL to reproject raster back into geographic coordinates, if needed
if geographic:
    deformation_raster_wgs84 = deformation_raster[:-3] + '_wgs84.tif'
    cmd = 'gdalwarp ' + deformation_raster + ' ' + deformation_raster_wgs84 + \
        ' -t_srs "+proj=longlat +ellps=WGS84"'
    call(cmd, shell=True)
             
