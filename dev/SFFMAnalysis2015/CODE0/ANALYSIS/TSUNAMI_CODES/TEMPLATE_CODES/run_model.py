"""Script for running a tsunami inundation scenario.

Source data such as elevation and boundary data is assumed to be available in
files specified by project.py

Run as:
> python runenbp.py
or  e.g. 
> python runenbp.py XX
or e.g. 
> mpirun -np 6 python runenbp.py XX > outfile.log 
where XX is an integer used to initiate a random_seed for the gallovic routine. 
If XX is left blank, then zero is used instead

Geoscience Australia, 2004-present
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
# Standard modules
import os
import time
import sys
import numpy
import scipy
import scipy.interpolate

# Related major packages
import anuga
from anuga.utilities import sww_merge

# Application specific imports
import project                 # Definition of file names and polygons
from anuga.tsunami_source import okada_tsunami

from anuga.parallel import distribute, myid, numprocs, finalize, barrier

from anuga import create_domain_from_regions


# Geometric pars
beach_slope=-1.0e-02
beach_end=5000.
shelf_start=0. # Shelf starts this far beyond the beach end
shelf_start_depth=0. # Shelf starts at this depth 
shelf_finish=200000. # Shelf ends this far beyond the beach end
shelf_bottom=2000. # Shelf drops this much
slope_finish=300000. # Slope ends this far beyond the beach
slope_bottom=6000. # Slope drops this much


# Build domain
if(myid==0):
    print 'Hello from processor ', myid
    # Make a directory to store outputs
    os.mkdir(project.output_run)
    print 'OUTPUT FOLDER IS: ', project.output_run
    anuga.copy_code_files(project.output_run,'run_model.py', 'project.py') # Write model files to folder
    
    #------------------------------------------------------------------------------
    # Create the triangular mesh and domain based on 
    # overall clipping polygon with a tagged
    # boundary and interior regions as defined in project.py
    #------------------------------------------------------------------------------
    anuga.create_mesh_from_regions(project.bounding_polygon,
                                   boundary_tags={'top': [0],
                                                  'ocean_east': [1],
                                                  'bottom': [2],
                                                  'onshore': [3]},
                               maximum_triangle_area=project.low_res_0,
                               filename=project.meshname,
                               interior_regions=[], #project.interior_regions,
                               breaklines=project.breakLines,
                               regionPtArea=project.regionPtAreas,
                               use_cache=False,
                               verbose=True)
    domain=anuga.create_domain_from_file(project.meshname)

    # Print some stats about mesh and domain
    print 'Number of triangles = ', len(domain)
    print 'The extent is ', domain.get_extent()
    print domain.statistics()
    
    domain.set_flow_algorithm('tsunami')
    domain.verbose=True

    #------------------------------------------------------------------------------
    # Setup parameters of computational domain
    #------------------------------------------------------------------------------
    domainName=project.name_stem + '_' + project.scenario
    domain.set_name(domainName) # Name of sww file
    domain.set_datadir(project.output_run)                       # Store sww output here
    
    #------------------------------------------------------------------------------
    # Setup initial conditions
    #------------------------------------------------------------------------------
    
    domain.set_quantity('friction', 0.02) 

    def topofun(x,y):
        # Idealised subduction zone topography
        L0=shelf_start #0.
        Z0=shelf_start_depth #0.
        L1=shelf_finish #200000.
        Z1=shelf_bottom #2000.
        L2=slope_finish #300000.
        Z2=slope_bottom #6000.
      
        S0=beach_slope
        S1=-Z1/L1
        S2=-(Z2-Z1)/(L2-L1) 
  
        # NOTE: min(x) = 0 in the coordinates passed to this function 
        xmin=beach_end #5000. # Offset allows positive topography near xmin
        topo=S0*(x-xmin)*(x-xmin < L0)+\
             (S1*(x-L0-xmin)+S0*L0)*(x-xmin>=L0)*(x-xmin <L1) +\
             (S2*(x-L1-xmin)+S1*L1)*(x-xmin >= L1)*(x-xmin < L2)+\
             (-Z2)*(x-xmin >= L2)
        return topo

    domain.set_quantity('elevation', topofun)
    
    domain.set_quantity('stage', project.tide)
    domain.set_quantity('xmomentum', 0.0)
    domain.set_quantity('ymomentum', 0.0)
else:
    domain= None
    print 'Hello from processor ', myid

domain=distribute(domain)

##------------------------------------------------------------------------------
## Setup information for earthquake scenario
##------------------------------------------------------------------------------
if project.scenario == 'earthquake':
  
    if(project.earthquake_type=='surface_perturb'):

        # Read xyz points describing the water surface perturbation.
        # Triangulate it (or use nearest neighbours for speed -- okay as we taper the edges of the deformation)
        sourceXYZ=scipy.genfromtxt(project.surface_perturbation_points,delimiter=',')
        # Linear interpolation is slow
        source_fun=scipy.interpolate.NearestNDInterpolator(sourceXYZ[:,0:2], sourceXYZ[:,2])

        def tsunami_source(x,y):
            # Need to map our xy points onto the along-strike /
            # along-surface-in-down-dip-direction coordinate system
            #
            # Approach: The near-surface along-strike boundary of the fault should align with the vertical line
            #           x=beach_end+slope_finish (where min(x) in the entire domain ==0)
            #           The centre of the fault should be on the line y=0
            newY=beach_end+slope_finish -x
            newX=(y+domain.geo_reference.yllcorner)
            outpts=source_fun(newX,newY)
            # Find nan values -- nan's are not equal to themeselves
            outpts_nan=scipy.nonzero(outpts!=outpts)
            outpts[outpts_nan]=0.
            return outpts
        
    else:
        err_mess= 'ERROR: earthquake type not recognized (check spelling)'
        raise sys.exit(err_mess) 

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------

# Try using the set-stage transmissive t zero n boundary condition here
def boundary_sealevel(t):
    return project.tide

Bt = anuga.shallow_water.boundaries.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(domain, boundary_sealevel)

# Boundary conditions for earthquake scenario
domain.set_boundary({'ocean_east': Bt,
                     'bottom': Bt,
                     'onshore': Bt,
                     'top': Bt})


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
t0 = time.time()

if project.scenario == 'earthquake':

    from anuga.operators.collect_max_stage_operator import Collect_max_stage_operator
    max_operator = Collect_max_stage_operator(domain)

    print ' Adding tsunami source to domain'
    print ' Adding to stage ...'
    domain.add_quantity('stage', tsunami_source)
    print ' Finished'

    print ' Adding to elevation ...'    
    domain.add_quantity('elevation', tsunami_source)
    print ' Finished'

    barrier()    
    if(myid==0):
        print 'Beginning main simulation'
   
    # Main evolve loop 
    for t in domain.evolve(yieldstep=project.yieldstep, 
                           finaltime=project.finaltime):
        print domain.timestepping_statistics()


if project.scenario == 'stationary':
    #print 'getting to the start of the loop'
    barrier()    
    for t in domain.evolve(yieldstep=project.yieldstep, 
                           finaltime=project.finaltime):
        print domain.timestepping_statistics()
        #print domain.boundary_statistics(tags='bottom')
        xx = domain.quantities['xmomentum'].centroid_values
        yy = domain.quantities['ymomentum'].centroid_values
        dd = domain.quantities['stage'].centroid_values - domain.quantities['elevation'].centroid_values
        dd = dd*(dd>0.)
        vv = (dd>1.0e-02)*( (xx/(dd+1.0e-06))**2 + (yy/(dd+1.0e-06))**2)**0.5
        print 'velocity max = ', vv.max()
#print 'That took %.2f seconds' %(time.time()-t0)

barrier()

# Get peak_stage outputs and save to file
c_v = max_operator.max_stage.centroid_values.flatten()
c_x = domain.centroid_coordinates[:,0].flatten()
c_y = domain.centroid_coordinates[:,1].flatten()
c_e = domain.quantities['elevation'].centroid_values.flatten()
out_max_data = numpy.vstack((c_x,c_y,c_v, c_e)).transpose()
out_max_file = project.output_run+"/" + 'max_stage_P' + str(myid)+'.txt'
numpy.savetxt(out_max_file, out_max_data,fmt='%.4f')

barrier()

if ((myid == 0) & (numprocs>1)):
    print 'Number of processors %g ' %numprocs
    print 'That took %.2f seconds' %(time.time()-t0)
    print 'Communication time %.2f seconds'%domain.communication_time
    print 'Reduction Communication time %.2f seconds'%domain.communication_reduce_time
    print 'Broadcast time %.2f seconds'%domain.communication_broadcast_time

    os.chdir(project.output_run)
    sww_merge.sww_merge_parallel(domainName, np=numprocs, verbose=True, delete_old=True)

finalize()
