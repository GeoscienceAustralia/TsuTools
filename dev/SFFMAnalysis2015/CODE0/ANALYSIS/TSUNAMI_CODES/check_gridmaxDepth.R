## Code to check the outputs of grid_maxDepth.py -- confirm it is working
##
## Idea: We look at a model run, read in the max-stage xyze,
##       then compute the depth, and compare with the wet-dry-line


# Location with max_stage*.txt files and wet_dry_line.txt
dataDir="/home/gareth/tmp/Nick_h_slip_paper/GARETH_TEST/TSUNAMI_MODELS/ANUGA_CODE/FOR_NCI/RUNS/S_1923009001_M7.95_Kanto-Japan_65/OceanInitial_1/MODEL_OUTPUTS/20140225_122248_run__earthquake_wsp_idealTsunami_missing_1000"

# Parameters from grid_maxDepth.py
dx=10.
dy=10.
width=10000.
length=10000.

# Read in max_stage*.txt files
ms=Sys.glob(paste0(dataDir,'/', 'max_stage*.txt'))
xyze=read.table(ms[1],sep=" ", header=FALSE)
for(i in 2:length(ms)) xyze=rbind(xyze, read.table(ms[i], sep=" ", header=FALSE))

# Compute max onshore depth
maxOnshoreDepth=(xyze[,3]-xyze[,4])*(xyze[,3]>xyze[,4]+1.0e-03)*(xyze[,4]>0.)

# Read wet-dry-line, made from contouring an image
wdLine=read.table(paste0(dataDir,'/','wet_dry_line.txt'))

# Convert from the 'image contour' that python gives, back to ANUGA coordinates
xmin=min(xyze[,1])
xmax=xmin+length
ymid=0.5*(min(xyze[,2]+max(xyze[,2])))
ymin=ymid - width/2.
ymax=ymid + width/2.

wdLine_2=wdLine
wdLine_2[,1]=(wdLine_2[,1]*dx+xmin)
wdLine_2[,2]=(wdLine_2[,2]*dy+ymin)

# Plot the transpose so I can see it better on my screen
# Lower y coordinate == more onshore
plot(wdLine_2[,2:1],t='l',asp=1,lwd=1,col='blue')
points(xyze[,2:1],col=c('red', 'black')[1+(maxOnshoreDepth>0.)],pch='.')
points(wdLine_2[,2:1],t='l',col='blue',lwd=2)

png('ExampleWDLineCheck.png',width=5,height=50,units='in',res=100)
plot(wdLine_2[,1:2],t='l',asp=1,lwd=1,col='blue')
points(xyze[,1:2],col=c('red', 'black')[1+(maxOnshoreDepth>0.)],pch='.')
points(wdLine_2[,1:2],t='l',col='blue',lwd=2)
dev.off()

