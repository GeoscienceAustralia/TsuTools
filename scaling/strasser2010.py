"""
Implements fault scaling relationship for subduction interface events
from Strasser et al. (2010).

Reference:
Strasser, F. O., Arango, M. C., & Bommer, J. J. (2010). Scaling of the source 
dimensions of interface and intraslab subduction-zone earthquakes with moment
magnitude. Seismological Research Letters, 81(6), 941-950.

Jonathan Griffin
Geoscience Australia
October 2014
"""

import numpy

def calculate_dimensions(magnitude):
    """
    Return fault dimensions (area, width, length) for input magnitude.
    :param magnitude
        Earthquake magnitude (Mw)
    """

    coefficients = {'length': [-2.477, 0.585],
                    'width': [-0.882, 0.351],
                    'area': [-3.476, 0.952]}
    length = numpy.power(10, (coefficients['length'][0] + \
                                  coefficients['length'][1]*magnitude)) 
    width = numpy.power(10, (coefficients['width'][0] + \
                                  coefficients['width'][1]*magnitude))
    area = numpy.power(10, (coefficients['area'][0] + \
                                  coefficients['area'][1]*magnitude))

    return length, width, area

def calculate_magnitude(dimension_type, dimension_value):
    """
    Return fault earthquake magnitude (Mw) for input fault dimensions 
    (area, width, length)
    :param dimension_type
        One of 'length', 'width', 'area'
    :param dimension_value
        value of the dimension in units of km (length, width) and km^2 (area)
    """

    coefficients = {'length': [4.868, 1.392],
                    'width': [4.410, 1.805],
                    'area': [4.441, 0.846]}
    magnitude = coefficients[dimension_type][0] + \
                    coefficients[dimension_type][1]*numpy.log10(dimension_value)

    return magnitude

