"""Creates and event by combining subfaults for a given
magnitude and centroid
Jonathan Griffin
February 2015
"""

import os, sys, errno
import numpy
from TsuTools.scaling import strasser2010
from TsuTools.utilities import plot_deformation
from subprocess import call

def nearest_point(point, points):
    """Function that finds the nearest point from a set
    of points, to a given point)
    point: 1D numpy array of len N
    points ND array of length M"""

    dist_2 = numpy.sum((points.transpose() - point)**2, axis=1)
    nearest_point_index = numpy.argmin(dist_2)
    return points.transpose()[nearest_point_index], nearest_point_index
                       


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
    data = numpy.genfromtxt(filename, skip_header = (num_headers+3))

    # Caculate dimensions
    length, width, area, sigma_length_plus, sigma_width_plus, \
        sigma_area_plus, simga_length_minus, sigma_width_minus, \
        sigma_area_minus = strasser2010.calculate_dimensions(magnitude)

    # Get required subfault information
    subfault_length = numpy.mean(data[:,7])
    subfault_width = numpy.mean(data[:,8])
    down_dip_number = int(max(data[:,10]))
    along_strike_number = int(max(data[:,9]))

    # Find nearest subfault to the index
    search_point = numpy.array([centroid_longitude, centroid_latitude])
    locations = numpy.array([data[:,0], data[:,1]])
    centroid_point, centroid_index = nearest_point(search_point, locations)
    centroid_strike_index = int(data[centroid_index,9])
    centroid_dip_index = int(data[centroid_index,10])
    
    fault_length = subfault_length*along_strike_number
    fault_width = subfault_width*down_dip_number

    # Find start index. First calculate half the length of the fault 
    # for the given magnitude, then how many subfaults this will be. 
    # Then count back from the centroid_strike_index this many subfault,
    # or until the edge of the fault is reached. Similarly for going up
    # dip based on the subfault width.
    num_subfaults_half_length = int(numpy.round((length/2)/subfault_length))
    num_subfaults_half_width = int(numpy.round((width/2)/subfault_width))
    print 'Area - Strasser (km^2)', area

    # If fault is not wide enough for number of subfaults, we can try 
    # adding extra subfaults in the along strike direction
    if 2*num_subfaults_half_width > down_dip_number:
        spare_dip_subfaults = 2*num_subfaults_half_width -\
            down_dip_number
        spare_area = area - (down_dip_number*subfault_width)* \
            (2*num_subfaults_half_length*subfault_length)
        spare_strike_subfaults = numpy.floor((spare_area/fault_width)/ \
                                                 subfault_length)
        # Adjust num_subfaults_half_length to extend along strike to match area
        num_subfaults_half_length += spare_strike_subfaults/2

    start_strike_index = int(max(0, centroid_strike_index - \
                                     num_subfaults_half_length))
    if start_strike_index == 0:
        remainder_subfaults = abs(centroid_strike_index - \
                                      num_subfaults_half_length)
    else:
        remainder_subfaults = 0
    end_strike_index = int(min(along_strike_number, centroid_strike_index + \
                                   num_subfaults_half_length + \
                                   remainder_subfaults))
    if end_strike_index == along_strike_number:
        remainder_subfaults_strike = abs(num_subfaults_half_length - \
                                      (along_strike_number - centroid_strike_index))
    else:
        remainder_subfaults_strike = 0
   
    start_dip_index = int(max(0, centroid_dip_index - num_subfaults_half_width))
    if start_dip_index == 0:
        remainder_subfaults = abs(centroid_dip_index - \
                                      num_subfaults_half_width)
    else:
        remainder_subfaults = 0
    end_dip_index = int(min(down_dip_number, centroid_dip_index + \
                                num_subfaults_half_width + \
                                remainder_subfaults))
    if end_dip_index == down_dip_number:
        remainder_subfaults_dip = abs(num_subfaults_half_width - \
                                          (down_dip_number - centroid_dip_index))
    else:
        remainder_subfaults_dip = 0
     
    # Try re-allocating spare subfaults back up strike and up dip
    # This should always work as earlier we check if there are 
    # enough subfaults down dip, but is not yet tested
    if remainder_subfaults_strike > 0:
        if start_strike_index > 0:
            start_strike_index = max(0, start_strike_index - remainder_subfaults_strike)
    if remainder_subfaults_dip > 0:
        if start_dip_index > 0:
            start_dip_index = max(0, start_dip_index - remainder_subfaults_dip)

    fault_area = (end_strike_index - start_strike_index)*subfault_length* \
        (end_dip_index - start_dip_index)*subfault_width

    if fault_area > sigma_area_plus or fault_area < sigma_area_minus:
        print 'Fault dimensions do not match scaling relationships for '\
            'given magnitude %.2f' % magnitude
        print 'Consider reducing the magnitude or the start index'
        sys.exit()

    # Calculate slip for actual fault area
    mu = 3e11 # dyne/cm assumed rigidity
    # Calculate moment
    moment_dynecm = numpy.power(10,((magnitude+10.73)*1.5))
    # Calculate slip
    slip = moment_dynecm/(mu*fault_area*1e10)
    slip = slip/100 # convert from cm to m
    print 'fault area', fault_area
    print 'slip', slip

    ################################################
    # Scale and add deformation for each subfault 
    ###############################################
    deformation_path = os.path.join(source_zone, 'deformation_data')
    start_grid_index = start_strike_index*down_dip_number + start_dip_index
    end_grid_index = end_strike_index*down_dip_number + end_dip_index
    start_grid_filename = os.path.join(deformation_path, '%s-%04d-1m.grd') % \
        (source_zone, start_grid_index)
    cmd = 'cp ' + start_grid_filename + ' tmp_sum.grd'
    call(cmd, shell=True)

    # Get indices of relevant subfaults
    index_list = []
    for i in range(start_strike_index, end_strike_index):
        for j in range(start_dip_index, end_dip_index):
            index_list.append(i*down_dip_number + j)
    print 'number of subfaults', len(index_list)
    for i in index_list:
        # Skip first subfaults, as already have it above
        if i == (start_strike_index*down_dip_number + start_dip_index):
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
        else: raise
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


        
    
