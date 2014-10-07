## The test elevation data was 'dropped' 100m with the following command

gdal_calc.py -A Res_500_elevation.tif --calc='A-100.0' --outfile 'Res_500_elevation_less_100.tif' 
