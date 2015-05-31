import os
import scipy
import scipy.interpolate
import glob
import matplotlib
matplotlib.use('Agg') # Avoid plot device errors on NCI
from matplotlib import pyplot

# Input parameters
#sourcedirs=glob.glob('EIGHT_FFI/*/*/*')
#sourcedirs = glob.glob('EIGHT_FFI/142782217353858_S_SC__NST_clip_SSD_none_RCS_TRUE_hurst_1_noise_stable/S_2011003011_M9_Tohoku-Oki--Japan_7/OceanInitial_9')
#sourcedirs=glob.glob('EIGHT_FFI2/*S_SCF*/*/*')
sourcedirs = glob.glob('BEST_MODEL_REMAINING_FFI/*/*/*')

dx = 10.
dy = 10.
width = 10000.
length = 10000.


for sourcedir in sourcedirs:
    basedir = sourcedir+'/MODEL_OUTPUTS/*/'

    # Read all the data
    max_stage_files = glob.glob(basedir+'max_stage_P*.txt')

    if len(max_stage_files)==0:
        continue

    for i, msFile in enumerate(max_stage_files):
        print 'Reading ', msFile
        tmpXYZ=scipy.genfromtxt(msFile)
        if i==0:
           fullXYZ=tmpXYZ
        else:
            fullXYZ=scipy.vstack((fullXYZ,tmpXYZ))

    # Compute depth, ensure positive (due to tsunami algorithm)
    myD=fullXYZ[:,2]-fullXYZ[:,3]
    myD=myD*(myD>0.)

    # Define extents of grid where we get outputs
    minX=min(fullXYZ[:,0])
    maxX=min(fullXYZ[:,0])+length
    minY=0.5*(min(fullXYZ[:,1])+max(fullXYZ[:,1]))-width/2.
    maxY=0.5*(min(fullXYZ[:,1])+max(fullXYZ[:,1]))+width/2.

    # Grid the data
    xinc=scipy.arange(minX,maxX,dx)
    yinc=scipy.arange(minY,maxY,dy)
    grid_x , grid_y = scipy.meshgrid(xinc,yinc)
    grid_depth= scipy.interpolate.griddata(fullXYZ[:,0:2], myD, 
                                           (grid_x, grid_y), method='nearest')
    if(False):
        # Might want this later
        grid_elev= scipy.interpolate.griddata(fullXYZ[:,0:2], fullXYZ[:,3], 
                                               (grid_x, grid_y), method='nearest')
        grid_stage= scipy.interpolate.griddata(fullXYZ[:,0:2], fullXYZ[:,2], 
                                               (grid_x, grid_y), method='nearest')

    # Extract the wet-dry line, assuming a single line. 
    wet_dry_contour=pyplot.contour(grid_depth>0., levels=[0.5])
    try:
        contour_path=wet_dry_contour.collections[0].get_paths()[0].vertices
        scipy.savetxt(os.path.dirname(msFile)+'/wet_dry_line.txt',contour_path)
    except:
        print 'Nontrivial contour at ', basedir
    if grid_depth.min() > 0.:
        print 'Full inundation: Contouring will fail so I am constructing the appropriate file'
        seq1 = scipy.linspace(0., 999., num=1000)
        # A contour at the very edge of the domain (full inundation)
        contour_path = scipy.vstack([seq1*0. + 0.5, seq1]).transpose()
        scipy.savetxt(os.path.dirname(msFile)+'/wet_dry_line.txt',contour_path)

