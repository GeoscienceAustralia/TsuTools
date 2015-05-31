# Quick R script to prepare input data
#

# Firstly, a function equivalent to the depth function in python is defined.
# This will allow us to compute the depth for smoothing, for the xyDef points
# Note that the function is defined relative to the orientation of the xyDef points
# (x==0 at base of trench, x increases towards land)

# Secondly, we read a file containing the dip/htop for the event
# We compute the x-shift Htop/(0.04+tan(dip)), which is the
#  amount we should move the deformation towards land, to make it occur on the
#  fault plane, and have the desired htop.

# Next, we compute the depth of the xyDef points

##################################################################
# Parameters from python
beach_slope=-1.0e-02
beach_end=5000.
shelf_start=0. # Shelf starts this far beyond the beach end
shelf_start_depth=0. # Shelf starts at this depth 
shelf_finish=200000. # Shelf ends this far beyond the beach end
shelf_bottom=2000. # Shelf drops this much
slope_finish=300000. # Slope ends this far beyond the beach
slope_bottom=6000. # Slope drops this much
    
trench_slope=(slope_bottom-shelf_bottom)/(slope_finish-shelf_finish)

depthFun<-function(x){
    # Here we code up the same depth function as in python,
    # Except viewed from the perspective of 'base of trench has x = 0'
    # x going towards shore = positive
    # This is the orientation of the deformation data
    #
    depth=x*0
    # Seaward of trench
    pts=which(x<0)
    depth[pts] = slope_bottom
    # On trench
    pts=which((x>=0) & (x<=(slope_finish-shelf_finish)))
    depth[pts]=slope_bottom-x[pts]*trench_slope
    # Landward of trench
    pts=which(x>(slope_finish-shelf_finish))
    depth[pts]=shelf_bottom - (x[pts] - (slope_finish-shelf_finish))*(shelf_start_depth-shelf_bottom)/(shelf_start-shelf_finish)

    # Don't allow negative depth
    depth=pmax(depth,0)

    return(depth)
}

####################################################################

# FIRST, we write the dip + htop to a simple file, and
# we will use the values to adjust the fault position
EqPar=read.table('INPUT_DATA/SrcMod_MetaData.csv',sep=",",header=T)

# Find _Number extension at the end of the source name
# This is at the end of the name of the directory above where we are
mywd=basename(dirname(getwd()))
myEventNum=as.numeric(strsplit(mywd,'_')[[1]][5])

# Extract dip + htop
dip=EqPar[myEventNum,]$DIP
htop=EqPar[myEventNum,]$Htop
# Compute the distance we need to push
xjust=htop/(trench_slope + tan(dip/180*pi))*1000
print(paste('Earthquake x offset: ', xjust))

write.table(c(dip, htop), 'INPUT_DATA/dip_htop.txt',sep=",",row.names=FALSE,col.names=FALSE)

##################################################

source('kajiura_filter.R')

INFILE=Sys.glob('INPUT_DATA/OceanInitial_*.xyz')
if(length(INFILE)!=1) stop('MORE THAN 1 INPUT DEFORMATION FILE FOUND')

xyDef=matrix(scan(INFILE,sep=","),ncol=3,byrow=T)

print('Adjusting the x coordinate of xyDef prior to smoothing')
xyDef[,1]=xyDef[,1]+xjust

print('Computing the static depth for each xyDef point')
depth=depthFun(xyDef[,1])

print('Smoothing')
# Apply the Kajiura filter
new_xyDef=tsunami_init_filter(xyDef, depth, grid_dx=1000, grid_dy=1000)

OUTFILE=gsub('.xyz','_SMOOTHED.xyz', INFILE)

write.table(new_xyDef, OUTFILE, sep=",", row.names=FALSE,col.names=FALSE)

print('FINISHED R PREPROCESSING')
