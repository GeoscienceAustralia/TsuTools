""" Common filenames and locations for topographic data, meshes and outputs.
    This file defines the parameters of the scenario you wish to run.
"""

import anuga
from time import localtime, strftime, gmtime
import numpy
import csv
import sys, os
from anuga.parallel import myid, numprocs, send, receive
    
#------------------------------------------------------------------------------
# Define scenario as either earthquake or stationary.
#------------------------------------------------------------------------------
#scenario = 'stationary'      # Still water
scenario = 'earthquake'       # Earthquake applied inside the domain

#------------------------------------------------------------------------------
# Define the region in which we have detailed polygons (to resolve inundation)
#------------------------------------------------------------------------------
detailed_region='idealTsunami'

#------------------------------------------------------------------------------------
# Define the earthquake forcing, for use if scenario='earthquake' (ignored otherwise)
#------------------------------------------------------------------------------------
#earthquake_type='synthetic_slip' # 'uniform_slip' or 'finite_inversion' or 'synthetic_slip'
#earthquake_type='uniform_slip'
#earthquake_type='finite_inversion'
earthquake_type='surface_perturb'

#---------------------------------------------------
# Tide = MSL
#---------------------------------------------------
tide = 0.0

#---------------------------------------------------
# Timestepping
#---------------------------------------------------
yieldstep=3600.*3.
finaltime=3600.*6.

#----------------------------------------------------------------------------------
#
# PARAMETERS BELOW HERE SHOULD NOT CHANGE FOR CASUAL USAGE
#
#----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
# Read from the command
# line 
#-----------------------------------------------------------------------------------
if(len(sys.argv)>1):
    random_seed_int=sys.argv[1] # Integer 
    eq_template=sys.argv[2] # character
else:
    random_seed_int=0
    eq_template='missing'

print 'project.random_seed_int', random_seed_int
print 'project.eq_template', eq_template


#------------------------------------------------------------------------------
# Get parameters for earthquake forcing
#------------------------------------------------------------------------------
if(scenario=='stationary'):
    # Useful flag for filenames
    run_nameflag="test"

if(earthquake_type=='surface_perturb'):
    run_nameflag='wsp'
    # Here we pass xyz points describing the ground/ocean perturbation
    surface_perturbation_points='INPUT_DATA/REPLACEWITHSED'


#------------------------------------------------------------------------------
# Create some filenames
#------------------------------------------------------------------------------
name_stem = scenario +'_'+run_nameflag+'_'+detailed_region+'_'+str(eq_template)+'_'+str(1000+random_seed_int)
print 'name_stem is '+ name_stem 


# Determine time for setting up local directories
if myid == 0:
    run_time=strftime('%Y%m%d_%H%M%S', localtime()) + '_run_'
    output_run='./MODEL_OUTPUTS/'+run_time+'_'+name_stem # Place to store outputs
    try:
        os.mkdir('./MODEL_OUTPUTS')
    except:
        pass

    for i in range(1, numprocs):
        send(output_run, i)
else:
    output_run = receive(0)


if(earthquake_type=='synthetic_slip'):
    # Location to save output synthetic earthqake files
    slip_output=output_run + '/' + synthetic_slip_dirname + '/' 

# Location to store the '.pts' file which is written prior to mesh generation
meshname = output_run + '/' + name_stem + '.msh'

#------------------------------------------------------------------------------
# Domain definitions
#------------------------------------------------------------------------------
# Use these very large zones as a boundary condition buffer.
bounding_polygon=anuga.read_polygon('INPUT_DATA/buffer_0.csv')
bounding_polygon_inner=anuga.read_polygon('INPUT_DATA/buffer_box.csv')
bounding_polygon_1 = anuga.read_polygon('INPUT_DATA/main_extent.csv')
# Site specific regions
if(detailed_region=='idealTsunami'):
    pass
else:
    err_mess='detailed_region not recognized'
    raise sys.exit(err_mess)

# Define resolutions (max area per triangle) for each polygon
#default_res = 400000*400000*0.5 # 400km triangles -- indended to be very very diffusive, to avoid boundary reflections
low_res_0=400000*400000*0.5 # 40km triangles -- intended to be very diffusive, avoid boundary reflections
low_res_1=40000*40000*0.5 # 40km triangles -- intended to be very diffusive, avoid boundary reflections
default_res_1 = 3000*3000*0.5 # Decent 'deep-water' resolution
interior_1_res=1000*1000*0.5 # Intermediate depth water resolution. 
interior_2_res=500.*500.*0.5 # Mod-shallow depth water resolution
interior_3_res=200.*200.*0.5 # Mod-shallow depth water resolution
interior_4_res=80.*80.*0.5 # shallow depth water resolution
interior_5_res=30.*30.*0.5 # shallow depth water resolution

## Define list of interior regions with associated resolutions
if(detailed_region=='none'):
    # Coarse domain for experimentation
    interior_regions = [
                        [bounding_polygon_1,default_res_1]]
elif(detailed_region=='idealTsunami'):
    poly1=anuga.read_polygon('INPUT_DATA/poly1.csv')
    poly2=anuga.read_polygon('INPUT_DATA/poly2.csv')
    poly3=anuga.read_polygon('INPUT_DATA/poly3.csv')
    poly4=anuga.read_polygon('INPUT_DATA/poly4.csv')
    poly5=anuga.read_polygon('INPUT_DATA/poly5.csv')
    #interior_regions = [
    #                    [bounding_polygon_inner,low_res_1],
    #                    [bounding_polygon_1,default_res_1],
    #                    [poly1, interior_1_res],
    #                    [poly2,interior_2_res],     
    #                    [poly3,interior_3_res],  
    #                    [poly4,interior_4_res], 
    #                    [poly5,interior_5_res]  
    #                    ]
    #
    breakLines=[bounding_polygon_inner,bounding_polygon_1, poly1, poly2,poly3, poly4, poly5]
    polyRes=[low_res_1, default_res_1, interior_1_res, interior_2_res, interior_3_res, interior_4_res, interior_5_res]
    regionPtAreas=[[0.5*(breakLines[i][1][0]+breakLines[i][2][0])-1.,\
                    0.5*(breakLines[i][1][1]+breakLines[i][2][1]),\
                    polyRes[i]] \
                    for i in range(len(breakLines))]
else:
    err_mess='detailed_region not recognized '
    raise sys.exit(err_mess)

