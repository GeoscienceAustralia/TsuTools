"""Creates and event by combining subfaults for a given
magnitude and centroid

Jonathan Griffin, Geoscience Australia, 2015
Gareth Davies, Geoscience Australia, 2015

"""

import os
import sys
import errno
import numpy
from TsuTools.scaling import strasser2010
from TsuTools.utilities import plot_deformation
from subprocess import call


def nearest_point(point, points):
    """Function that finds the nearest point from a set
    of points, to a given point)
    point: 1D numpy array of len N (N = 2 for 2 dimensional data)
    points: numpy array of shape (N, M) where M is the number of points

    return: The entry in points nearest to point, and its index"""

    # Note this uses euculidean distances -- so beware possible inaccuracy
    # using it on geographic coordinates at high latitudes. (Not sure how
    # extreme the situation has to be for it to matter -- does it ever?)
    dist_2 = numpy.sum((points.transpose() - point) ** 2, axis=1)
    nearest_point_index = numpy.argmin(dist_2)
    return points.transpose()[nearest_point_index], nearest_point_index


def distance_points_sphere(points, point, radius=6371.0):
    """Find the distance between each row of points with point.
        All coordinates are assumed to be lon/lat on a sphere with given radius

        points: numpy array of size (N, 2) with points[:,0] having longitude, and points[:,1] latitude
                Beware: This array orientation is different to the other routines here
        point:  numpy array of size (2,) containing [lon, lat]
        radius: Spherical radius (default is radius of the earth)

        return: The distances
    """
    # Convert coordinates to radians
    deg2rad = 2.0 * numpy.pi / 360.0
    lon1 = point[0] * deg2rad
    lat1 = point[1] * deg2rad
    lon2 = points[:, 0] * deg2rad
    lat2 = points[:, 1] * deg2rad

    # Haversine formula for central angle dtheta
    dlat = lat1 - lat2
    dlon = lon1 - lon2
    tmp = (
        numpy.sin(
            dlat /
            2.0) ** 2 +
        numpy.cos(lat1) *
        numpy.cos(lat2) *
        numpy.sin(
            dlon /
            2.0) ** 2) ** 0.5
    dtheta = 2.0 * numpy.arcsin(tmp)

    distance = radius * dtheta

    return(distance)


def nearest_point_sphere(point, points, radius=6371.0):
    """Function that finds the nearest point from a set
    of points, to a given point). All points assumed to be lon/lat on a sphere
    with prescribed radius.
    point: numpy array of length 2 containing lon, lat
    points: 2D numpy array of size (2,M) where M is the number of points.
            points[0,:] has lon, points[1,:] has lat
    radius: Sphere radius (default = radius of earth in km)

    return: The nearest point, and its index in points"""

    distance = distance_points_sphere(points.transpose(), point, radius=radius)

    nearest_point_index = numpy.argmin(distance)

    return points.transpose()[nearest_point_index], nearest_point_index


def get_subfaults_in_envelope(
        centroid_strike_number,
        centroid_dip_number,
        strike_numbers,
        dip_numbers,
        strike_envelope,
        dip_envelope):
    """Suppose subfaults are all defined by their along-strike and along-dip
    numbers (only).  These numbers correspond to cells in a logically
    rectangular array -- but some cells may be missing.  Given the 'centroid'
    strike and dip numbers, we find all subfaults within the 'envelope' defined by
    'strike_evelope, dip_envelope'. Subfaults with a strike and dip number that
    differ from the centroid numbers by an amount in the range of strike_envelope,
    dip_envelope are 'in the envelope'.

    centroid_strike_number:
    centroid_dip_number:
    strike_numbers: Array with strike number for all subfaults
    dip_numbers: Array with dip number for all subfaults
    strike_envelope: Array of 2 integers, giving negative/positive strike
        number difference inside the envelope
    dip_envelope: Array of 2 integers, giving negative/positive dip number
        difference inside the envelope

    return: The indices of all subfaults inside the envelope
    """

    assert len(strike_envelope) == 2, 'Strike envelope must have length = 2'
    assert len(dip_envelope) == 2, 'Dip envelope must have length = 2'

    assert strike_envelope[
        0] <= 0, 'First entry of strike envelope must be <= 0'
    assert dip_envelope[0] <= 0, 'First entry of dip envelope must be <= 0'
    assert strike_envelope[
        1] >= 0, 'Second entry of strike envelope must be >= 0'
    assert dip_envelope[1] >= 0, 'Second entry of dip envelope must be >= 0'

    dstrike = strike_numbers - centroid_strike_number
    ddip = dip_numbers - centroid_dip_number

    strike_indices_in_envelope = (dstrike >= strike_envelope[0]) *\
        (dstrike <= strike_envelope[1])
    dip_indices_in_envelope = (ddip >= dip_envelope[0]) *\
        (ddip <= dip_envelope[1])

    inside_envelope = strike_indices_in_envelope * dip_indices_in_envelope

    return inside_envelope.nonzero()[0]


def find_subfaults_to_include(
        along_strike_numbers,
        down_dip_numbers,
        centroid_index,
        desired_subfaults_length,
        desired_subfaults_width,
        prefer_down_dip_expansion=0,
        prefer_along_strike_expansion=0):
    """Suppose subfaults are arranged in a logically rectangular array (but
    with missing subfaults allowed) defined by the along_strike_numbers and
    along_dip_numbers. We want to make a rupture with approximate centroid
    at centroid_index, covering a set of subfaults with dimensions [on the
    logical grid] length = desired_subfaults_length, width = desired_subfaults_width.

    To do this we find all subfaults associated with increasingly wide strike
    and dip envelopes around the centroid subfault (see
    get_subfaults_in_envelope). We keep expanding the strike and dip envelopes,
    as long as they continue to reduce some measure of the difference between
    the desired shape of the rupture and its area.

    The way this is coded (function minimization) might help if we want to
    try different criteria for choosing the subfaults to include.

    The result can be affected by whether we expand first in the along strike
    (vs reverse strike) directions, and similarly for the dip. This can be
    controlled by the provided parameters

    along_strike_numbers: Integer array. Along-strike coordinates of
        subfaults in the logically rectangular array
    down_dip_numbers: Integer array (same length as along_strike_numbers). Down
        dip coordinates of subfaults in the logically rectangular array
    centroid_index: Index of desired centroid in 'along_strike_numbers' and
        'down_dip_numbers'.
    desired_subfaults_length: The rupture should span this range of along-strike
        subfaults in the logically rectangular layout
    desired_subfaults_width: The rupture should span this range of down-dip
        subfaults in the logically rectangular layout
    prefer_down_dip_expansion: If 0, try to extend up dip first, then down dip.
        If 1, try to extend down dip first, then up dip
    prefer_along_strike_expansion: If 0, try to extend in the 'reverse' strike
        direction first, then along the strike. If 1, try to extend along the
        strike first, then in the reverse strike direction.

    return: The indices of the subfaults to include in the event
    """

    assert (prefer_down_dip_expansion in [0, 1])
    assert (prefer_along_strike_expansion in [0, 1])

    # First make a function defining the discrepancy between a proposal
    # of subfaults to include, and the desired rupture properties
    def function_to_minimize(subfaults_included):
        """This function measures the discrepency between the width/length of
           a hypothetical rupture [defined by a set of subfaults] and
           the ideal values defined above

           It is used to help numerically search for the best set of subfaults

        """
        included_strike_numbers = along_strike_numbers[subfaults_included]
        included_dip_numbers = down_dip_numbers[subfaults_included]
        included_subfaults_length = included_strike_numbers.max() -\
            included_strike_numbers.min() + 1
        included_subfaults_width = included_dip_numbers.max() -\
            included_dip_numbers.min() + 1

        desired_subfaults_area = desired_subfaults_length * \
            desired_subfaults_width

        residual = abs(included_subfaults_length - desired_subfaults_length) +\
            abs(included_subfaults_width - desired_subfaults_width) +\
            abs(len(subfaults_included) - desired_subfaults_area) * 0

        return residual

    # Starting from a rupture defined by only the centroid index, grow
    # the rupture in all directions to reduce 'function_to_minimize'
    strike_envelope = numpy.array([0, 0])
    dip_envelope = numpy.array([0, 0])
    residual_last = 1e+30  # Larger than any value we will every compute (?)
    residual_start = function_to_minimize(numpy.array([centroid_index]))
    residual_now = residual_start

    centroid_strike_number = along_strike_numbers[centroid_index]
    centroid_dip_number = down_dip_numbers[centroid_index]

    # FIXME: This is a  crude search algorithm. There must be better methods?
    while residual_last > residual_now:
        residual_last = 1 * residual_now

        # Change the strike envelope
        search_order = [prefer_along_strike_expansion,
                        1 - prefer_along_strike_expansion]
        for i in search_order:
            strike_envelope_temp = 1 * strike_envelope
            strike_envelope_temp[i] += int(2.0 * (i - 0.5))  # -1 or +1
            subfaults_temp = get_subfaults_in_envelope(
                centroid_strike_number,
                centroid_dip_number,
                along_strike_numbers,
                down_dip_numbers,
                strike_envelope_temp,
                dip_envelope)
            new_resid = function_to_minimize(subfaults_temp)
            if new_resid < residual_now:
                strike_envelope = strike_envelope_temp
                residual_now = new_resid

        # Change the dip envelope
        search_order = [
            prefer_down_dip_expansion,
            1 -
            prefer_down_dip_expansion]
        for i in search_order:
            dip_envelope_temp = 1 * dip_envelope
            dip_envelope_temp[i] += int(2.0 * (i - 0.5))  # -1 or +1
            subfaults_temp = get_subfaults_in_envelope(
                centroid_strike_number,
                centroid_dip_number,
                along_strike_numbers,
                down_dip_numbers,
                strike_envelope,
                dip_envelope_temp)
            new_resid = function_to_minimize(subfaults_temp)
            if new_resid < residual_now:
                dip_envelope = dip_envelope_temp
                residual_now = new_resid

        subfaults_included = get_subfaults_in_envelope(
            centroid_strike_number,
            centroid_dip_number,
            along_strike_numbers,
            down_dip_numbers,
            strike_envelope,
            dip_envelope)

    return subfaults_included


def create_event_new(source_zone, magnitude,
                     centroid_longitude, centroid_latitude,
                     prefer_down_dip_expansion=0,
                     prefer_along_strike_expansion=0):
    """Function that converts i_invall file to Okada format

    """

    # Read data
    filename = os.path.join(source_zone, 'i_invall-%s' % source_zone)
    f_in = open(filename, 'r')
    fault_name = f_in.readline()
    msg = 'Source zone %s does not match name within i_invall '\
        'file %s' % (source_zone, fault_name)

    num_headers = int(f_in.readline())
    f_in.close()
    data = numpy.genfromtxt(filename, skip_header=(num_headers + 3))

    # Caculate dimensions from Strasser 2010 scaling relations theory
    length, width, area, sigma_length_plus, sigma_width_plus, \
        sigma_area_plus, sigma_length_minus, sigma_width_minus, \
        sigma_area_minus = strasser2010.calculate_dimensions(magnitude)

    print 'Strasser scaling relation results for Mw: ', magnitude
    print 'Area (km^2)', area
    print 'length (km): ', length
    print 'width (km): ', width
    print 'length x width (km^2) ( ...not exactly equal to Area... )',\
          length * width

    # Get required subfault information
    subfault_length = numpy.mean(data[:, 7])
    subfault_width = numpy.mean(data[:, 8])
    down_dip_number = int(max(data[:, 10]))
    along_strike_number = int(max(data[:, 9]))

    subfault_lengths = data[:, 7]
    subfault_widths = data[:, 8]
    along_strike_numbers = data[:, 9]
    down_dip_numbers = data[:, 10]

    # Find nearest subfault to the index
    search_point = numpy.array([centroid_longitude, centroid_latitude])
    locations = data[:, 0:2]
    centroid_point, centroid_index = nearest_point_sphere(search_point,
                                                          locations)
    centroid_strike_number = int(data[centroid_index, 9])
    centroid_dip_number = int(data[centroid_index, 10])

    fault_length = subfault_length * along_strike_number
    fault_width = subfault_width * down_dip_number

    # This defines approximately how long we would like the rupture to be
    num_subfaults_length = int(numpy.round((length * 1.0) / subfault_length))
    num_subfaults_width = int(numpy.round((width * 1.0) / subfault_width))

    # Here, find subfaults in the envelope, and get a measure of length, width
    # and area. Grow the envelope to find the 'best' agreement with these by
    # some measure.
    subfaults_included = find_subfaults_to_include(
        along_strike_numbers,
        down_dip_numbers,
        centroid_index,
        desired_subfaults_length=num_subfaults_length,
        desired_subfaults_width=num_subfaults_width,
        prefer_down_dip_expansion=prefer_down_dip_expansion,
        prefer_along_strike_expansion=prefer_along_strike_expansion)

    #
    # At this stage subfaults_included should be optimal
    #

    fault_area = (subfault_lengths[subfaults_included] *
                  subfault_widths[subfaults_included]).sum()

    if (fault_area > sigma_area_plus) or (fault_area < sigma_area_minus):
        print 'Fault dimensions do not match scaling relationships for '\
            'given magnitude %.2f' % magnitude
        print 'Consider reducing the magnitude or the start index'
        sys.exit()

    # Calculate slip for actual fault area
    mu = 3.0e11  # dyne/cm assumed rigidity
    # Calculate moment
    moment_dynecm = numpy.power(10.0, ((magnitude + 10.73) * 1.5))
    # Calculate slip
    slip = moment_dynecm / (mu * fault_area * 1.0e10)
    slip = slip / 100.0  # convert from cm to m
    print 'fault area', fault_area
    print 'slip', slip

    ################################################
    # Scale and add deformation for each subfault
    ###############################################
    deformation_path = os.path.join(source_zone, 'deformation_data')
    start_grid_index = subfaults_included[0]
    start_grid_filename = os.path.join(deformation_path, '%s-%04d-1m.grd') % \
        (source_zone, start_grid_index)
    cmd = 'cp ' + start_grid_filename + ' tmp_sum.grd'
    call(cmd, shell=True)

    # Get indices of relevant subfaults
    print 'number of subfaults', len(subfaults_included)
    for i in subfaults_included[1:]:
        grd_file = os.path.join(deformation_path, '%s-%04d-1m.grd') % \
            (source_zone, i)
        cmd = 'grdmath tmp_sum.grd ' + grd_file + ' ADD = tmp_sum.grd'
        print cmd
        call(cmd, shell=True)
    # Now scale grids by slip and put in a directory for events
    event_dir = source_zone + '/events'
    try:
        os.makedirs(event_dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(event_dir):
            pass
        else:
            raise
    print magnitude, centroid_longitude, centroid_latitude
    event_grd_name = '%s_Mw%.2f_%.3f_%.3f.grd' % (source_zone, magnitude,
                                                  centroid_longitude,
                                                  centroid_latitude)
    event_grd_path = os.path.join(event_dir, event_grd_name)
    cmd = 'grdmath tmp_sum.grd ' + str(slip) + ' MUL = ' + event_grd_path
    print cmd
    call(cmd, shell=True)
    # Convert to ESRI ascii raster
    event_ascii_path = event_grd_path[:-3] + 'asc'
    cmd = 'grdreformat ' + event_grd_path + ' ' + event_ascii_path + '=ef -V'
    print cmd
    call(cmd, shell=True)

    ##########################################
    # Now plot the deformation
    #########################################
    plot_deformation.plot_deformation(event_grd_path)


def create_event(source_zone, magnitude,
                 centroid_longitude, centroid_latitude):
    """Function that converts i_invall file to okada format
    """
    # Read data
    filename = os.path.join(source_zone, 'i_invall-%s' % source_zone)
    f_in = open(filename, 'r')
    fault_name = f_in.readline()
    msg = 'Source zone %s does not match name within i_invall '\
        'file %s' % (source_zone, fault_name)

    num_headers = int(f_in.readline())
    f_in.close()
    data = numpy.genfromtxt(filename, skip_header=(num_headers + 3))

    # Caculate dimensions
    length, width, area, sigma_length_plus, sigma_width_plus, \
        sigma_area_plus, simga_length_minus, sigma_width_minus, \
        sigma_area_minus = strasser2010.calculate_dimensions(magnitude)

    # Get required subfault information
    subfault_length = numpy.mean(data[:, 7])
    subfault_width = numpy.mean(data[:, 8])
    down_dip_number = int(max(data[:, 10]))
    along_strike_number = int(max(data[:, 9]))

    # Find nearest subfault to the index
    search_point = numpy.array([centroid_longitude, centroid_latitude])
    locations = numpy.array([data[:, 0], data[:, 1]])
    centroid_point, centroid_index = nearest_point(search_point, locations)
    centroid_strike_index = int(data[centroid_index, 9])
    centroid_dip_index = int(data[centroid_index, 10])

    fault_length = subfault_length * along_strike_number
    fault_width = subfault_width * down_dip_number

    # Find start index. First calculate half the length of the fault
    # for the given magnitude, then how many subfaults this will be.
    # Then count back from the centroid_strike_index this many subfault,
    # or until the edge of the fault is reached. Similarly for going up
    # dip based on the subfault width.
    num_subfaults_half_length = int(
        numpy.round(
            (length / 2) / subfault_length))
    num_subfaults_half_width = int(numpy.round((width / 2) / subfault_width))
    print 'Area - Strasser (km^2)', area

    # If fault is not wide enough for number of subfaults, we can try
    # adding extra subfaults in the along strike direction
    if 2 * num_subfaults_half_width > down_dip_number:
        spare_dip_subfaults = 2 * num_subfaults_half_width -\
            down_dip_number
        spare_area = area - (down_dip_number * subfault_width) * \
            (2 * num_subfaults_half_length * subfault_length)
        spare_strike_subfaults = numpy.floor((spare_area / fault_width) /
                                             subfault_length)
        # Adjust num_subfaults_half_length to extend along strike to match area
        num_subfaults_half_length += spare_strike_subfaults / 2

    start_strike_index = int(max(0, centroid_strike_index -
                                 num_subfaults_half_length))
    if start_strike_index == 0:
        remainder_subfaults = abs(centroid_strike_index -
                                  num_subfaults_half_length)
    else:
        remainder_subfaults = 0
    end_strike_index = int(min(along_strike_number, centroid_strike_index +
                               num_subfaults_half_length +
                               remainder_subfaults))
    if end_strike_index == along_strike_number:
        remainder_subfaults_strike = abs(
            num_subfaults_half_length - (along_strike_number - centroid_strike_index))
    else:
        remainder_subfaults_strike = 0

    start_dip_index = int(
        max(0, centroid_dip_index - num_subfaults_half_width))
    if start_dip_index == 0:
        remainder_subfaults = abs(centroid_dip_index -
                                  num_subfaults_half_width)
    else:
        remainder_subfaults = 0
    end_dip_index = int(min(down_dip_number, centroid_dip_index +
                            num_subfaults_half_width +
                            remainder_subfaults))
    if end_dip_index == down_dip_number:
        remainder_subfaults_dip = abs(num_subfaults_half_width -
                                      (down_dip_number - centroid_dip_index))
    else:
        remainder_subfaults_dip = 0

    # Try re-allocating spare subfaults back up strike and up dip
    # This should always work as earlier we check if there are
    # enough subfaults down dip, but is not yet tested
    if remainder_subfaults_strike > 0:
        if start_strike_index > 0:
            start_strike_index = max(
                0,
                start_strike_index -
                remainder_subfaults_strike)
    if remainder_subfaults_dip > 0:
        if start_dip_index > 0:
            start_dip_index = max(0, start_dip_index - remainder_subfaults_dip)

    fault_area = (end_strike_index - start_strike_index) * subfault_length * \
        (end_dip_index - start_dip_index) * subfault_width

    if fault_area > sigma_area_plus or fault_area < sigma_area_minus:
        print 'Fault dimensions do not match scaling relationships for '\
            'given magnitude %.2f' % magnitude
        print 'Consider reducing the magnitude or the start index'
        sys.exit()

    # Calculate slip for actual fault area
    mu = 3e11  # dyne/cm assumed rigidity
    # Calculate moment
    moment_dynecm = numpy.power(10, ((magnitude + 10.73) * 1.5))
    # Calculate slip
    slip = moment_dynecm / (mu * fault_area * 1e10)
    slip = slip / 100  # convert from cm to m
    print 'fault area', fault_area
    print 'slip', slip

    ################################################
    # Scale and add deformation for each subfault
    ###############################################
    deformation_path = os.path.join(source_zone, 'deformation_data')
    start_grid_index = start_strike_index * down_dip_number + start_dip_index
    end_grid_index = end_strike_index * down_dip_number + end_dip_index
    start_grid_filename = os.path.join(deformation_path, '%s-%04d-1m.grd') % \
        (source_zone, start_grid_index)
    cmd = 'cp ' + start_grid_filename + ' tmp_sum.grd'
    call(cmd, shell=True)

    # Get indices of relevant subfaults
    index_list = []
    for i in range(start_strike_index, end_strike_index):
        for j in range(start_dip_index, end_dip_index):
            index_list.append(i * down_dip_number + j)
    print 'number of subfaults', len(index_list)
    for i in index_list:
        # Skip first subfaults, as already have it above
        if i == (start_strike_index * down_dip_number + start_dip_index):
            continue
        grd_file = os.path.join(deformation_path, '%s-%04d-1m.grd') % \
            (source_zone, i)
        cmd = 'grdmath tmp_sum.grd ' + grd_file + ' ADD = tmp_sum.grd'
        print cmd
        call(cmd, shell=True)
    # Now scale grids by slip and put in a directory for events
    event_dir = source_zone + '/events'
    try:
        os.makedirs(event_dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(event_dir):
            pass
        else:
            raise
    print magnitude, centroid_longitude, centroid_latitude
    event_grd_name = '%s_Mw%.2f_%.3f_%.3f.grd' % (source_zone, magnitude,
                                                  centroid_longitude,
                                                  centroid_latitude)
    event_grd_path = os.path.join(event_dir, event_grd_name)
    cmd = 'grdmath tmp_sum.grd ' + str(slip) + ' MUL = ' + event_grd_path
    print cmd
    call(cmd, shell=True)
    # Convert to ESRI ascii raster
    event_ascii_path = event_grd_path[:-3] + 'asc'
    cmd = 'grdreformat ' + event_grd_path + ' ' + event_ascii_path + '=ef -V'
    print cmd
    call(cmd, shell=True)

    ##########################################
    # Now plot the deformation
    #########################################
    plot_deformation.plot_deformation(event_grd_path)
#    cmd = 'python plot_deformation.py %s' % event_grd_path
#    print cmd
#    call(cmd, shell=True)


if __name__ == "__main__":

    def usage():
        print 'python generate_event_deformation.py <source_zone> <Mw> '\
            '<centroid_longitude> <centroid_latitude>'

    try:
        filename = sys.argv[1]
        magnitude = float(sys.argv[2])
    except:
        usage()
        sys.exit()
    try:
        centroid_longitude = float(sys.argv[3])
        centroid_latitude = float(sys.argv[4])
    except:
        #start_index = 0
        usage()
        sys.exit()

    create_event(filename, magnitude, centroid_longitude, centroid_latitude)
