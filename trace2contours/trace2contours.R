## Convert a fault trace shapefile to a set of contours representing the fault
## plane

suppressPackageStartupMessages(library(methods))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rgdal))
suppressPackageStartupMessages(library(sp))
suppressPackageStartupMessages(library(geosphere))
suppressPackageStartupMessages(library(raster))

#' Extend a fault trace to fault-plane contours
#'  
#' Given a fault trace known to be dipping to the right of the line
#' orientation, with dip provided in the attribute table, and a known
#' maximum depth, extend 'normals' from the fault trace down the dip
#' to the max depth. We can then use this to make a polygon / interpolate
#' from the trace
#'
#' @param fault_trace SpatialLines of the fault trace, with appropriate
#'                      attribute table. 
#' @param maxdepth max depth of the slab, in KM
#' @param mindepth min depth of the slab, in KM
#' @param contour_maxdepth max output contour depth (should be < maxdepth)
#' @param contour_mindepth min output contour depth (should be > mindepth)
#' @param contour_interval (approximate) interval for output contours
#' @param in_0_360 Is the fault trace longitude within (0-360)
#' @param verbose 
#' @param plot_delay If verbose, keep plot on screen for 'plot_delay' seconds
#' @return contours from slab depths

extend_trace_to_depth_contours<-function(fault_trace, 
                                         maxdepth, 
                                         mindepth, 
                                         contour_maxdepth,
                                         contour_mindepth,
                                         contour_interval=2,
                                         in_0_360=TRUE,
                                         verbose=TRUE,
                                         plot_delay=0.0
                                         ){

    if(maxdepth > 120){
        stop('maxdepth>120, but I have assumed maxdepth is in km')
    }

    if(maxdepth < mindepth) stop('maxdepth < mindepth !')

    if(contour_mindepth <= mindepth){
         stop(paste0('contour_mindepth must be > fault_trace_depth'))
    }

    # Get points along the trace, and remove repeated points
    fault_trace_points = as(fault_trace, 'SpatialPointsDataFrame')
    keep_points=c(seq(1, length(fault_trace_points), 2), 
                  length(fault_trace_points))
    fault_trace_points = coordinates(fault_trace_points)[keep_points,]
    

    # Compute bearing of down dip direction, and find point at the maximal
    # depth along this bearing
    store_seg_destpoints=c()
    for( i in 1:length(fault_trace) ){
        segment_pts = as(fault_trace[i, ], 'SpatialPointsDataFrame')

        # Check that there are just 2 segment points
        if(length(segment_pts)!= 2){
            stop('Fault segment not defined by a single line')
        }

        fault_bearing = bearingRhumb(segment_pts[1,], segment_pts[2,])
        fault_downdip_bearing = (fault_bearing+90)%%360 # Normal, in [0-360]

        depth_m = (maxdepth - mindepth)*1000 # Depth in metres
        max_distance = depth_m/tan(segment_pts[1, ]@data$Dip/180*pi)

        seg_destpoint1 = destPoint(segment_pts[1,], fault_downdip_bearing, 
            d=max_distance)         
        if(in_0_360 & (seg_destpoint1[1] < 0) ){
            seg_destpoint1[1] = 360 + seg_destpoint1[1]
        }

        seg_destpoint2 = destPoint(segment_pts[2,], fault_downdip_bearing, 
            d=max_distance)         
        if( in_0_360 & (seg_destpoint2[1] < 0) ){
            seg_destpoint2[1] = 360 + seg_destpoint2[1]
        }

        store_seg_destpoints = rbind(store_seg_destpoints, seg_destpoint1)
        store_seg_destpoints = rbind(store_seg_destpoints, seg_destpoint2)

    }

    # We now have 2 seg_destpoints for each fault trace point, except at the
    # end points. Let's use the arithmetic average to remove the repeated
    # points, and thus determine the fault boundary
    ll = length(fault_trace) + 1
    thinned_seg_destpoints = matrix(NA, nrow=ll, ncol=2)
    for(i in 1:ll){
    
        if(i == 1){
            # Special case -- first point 
            thinned_seg_destpoints[i,] = store_seg_destpoints[1,]

        }else if(i == ll){
            # Special case -- last point
            # Get the last seg-destpoint. There are 2*(length(fault_trace))
            # in total
            thinned_seg_destpoints[i,] = store_seg_destpoints[2*i-2,]
        }else{
            # Simple arithmetic average of the normals extended from the point
            thinned_seg_destpoints[i,] = midPoint(store_seg_destpoints[2*i-2,],
                                                  store_seg_destpoints[2*i-1,])

            if(in_0_360 & (thinned_seg_destpoints[i,1] < 0)){
                thinned_seg_destpoints[i,1] = 360 + thinned_seg_destpoints[i,1]
            }

        }
    }
    
    # Make a polygon surrounding the fault in areas where we have depths
    fault_seg_poly=rbind(fault_trace_points, thinned_seg_destpoints[ll:1,])
    fault_pt_depths=c( rep(mindepth,ll), 
                       rep(maxdepth,ll))
    fault_seg_poly_closed=rbind(fault_seg_poly, fault_seg_poly[1,]) # Close ring
    fault_seg_poly_sp=
        SpatialPolygons(list(Polygons(
                             list(Polygon(fault_seg_poly_closed, hole=FALSE))
                                ,ID='0')),
                        proj4string=CRS('+init=epsg:4326'))

    # Now, we want to interpolate these depths onto a raster, then contour
    # 
    # Inefficient good-enough approach: Let's generate loads of long,lat,depth
    # points, then just do nearest neighbour or something
    if(verbose){
        x11() # Open graphics device for plotting
        plot(fault_seg_poly_sp, col='green', axes=T, asp=1)
        points(fault_seg_poly, col=as.factor(fault_pt_depths))
    }
    store_pts = c()
    box_horiz_n = 100
    box_down_n = 100
    for(i in 1:length(fault_trace)){
        box_start = fault_seg_poly[i:(i+1),]
        box_end = fault_seg_poly[(2*ll-(i-1)):(2*ll-1-(i-1)),]

        # Plot it
        mybox = rbind(box_start, box_end[2:1,])
        if (verbose) points(mybox, t='l')
        
        for(j in 1:box_horiz_n){

            step = (j - 1)/(box_horiz_n - 1)
            bp1 = box_start[1,] + (box_start[2,] - box_start[1,])*step
            bp2 = box_end[1,] + (box_end[2,] - box_end[1,])*step

            local_pts=cbind(
                gcIntermediate(bp1, bp2, n=box_down_n-2, addStartEnd=TRUE),
                seq(mindepth, maxdepth, len=box_down_n) )

            if(in_0_360 & any(local_pts[,1] < 0)){
                local_pts[,1] = local_pts[,1]*(local_pts[,1] >= 0) + 
                    (360 + local_pts[,1])*(local_pts[,1] < 0)
            }
            store_pts=rbind(store_pts, local_pts)
        }
    }

    # Push onto raster, in a crude way
    myrast = raster(extent(fault_seg_poly_sp)) #,nrow=50,ncol=50)
    res(myrast) = c(1, 1)/10 # 1/10 of a degree
    proj4string(myrast) = CRS('+init=epsg:4326')
    myrast = setValues(myrast, NA)

    store_pts_coords = SpatialPoints(store_pts[,1:2], 
        proj4string=CRS('+init=epsg:4326'))

    store_pts_indices = extract(myrast, store_pts_coords, cellnumbers=T)[,1]

    rast_vals = aggregate(store_pts[,3], list(store_pts_indices), FUN=mean,
        na.rm=T)

    myrast[rast_vals[,1]] = rast_vals[,2]
    myrast=extend(myrast, 10) # Add a buffer so that contouring works ok

    # Now, let's ensure the trace is burned in with the mindepth
    trace_pt_coords = extract(myrast, 
                              spsample(fault_trace,n=1000,type='regular'),
                              small=T, 
                              cellnumbers=T)
    trace_pt_coords = trace_pt_coords[,1]
    myrast[trace_pt_coords] = mindepth
    
    if(verbose){
        plot(myrast, col=rev(terrain.colors(255)), alpha=0.85, add=T, 
            zlim=c(contour_mindepth, contour_maxdepth))
    }

    # Make contour
    #contourLevels = c(mindepth + 0.001, 
    #                  seq(mindepth + contour_interval, 
    #                      maxdepth - contour_interval,
    #                      by = contour_interval), 
    #                  maxdepth - 0.001)
    contourLevels = seq(contour_mindepth, contour_maxdepth, 
        len=ceiling((contour_maxdepth-contour_mindepth)/contour_interval))

    mycontours = rasterToContour(myrast, level = contourLevels)

    plot(mycontours, add=T)

    output = list(fault_seg_poly_sp, myrast, mycontours)

    if(verbose) Sys.sleep(5)

    return(output)

}

#' Ensure contours are oriented with the dip to the right, like the trace
#'
#' @param contours fault trace contours with potentially incorrect orientation
#' @param fault_trace fault trace with correct orientation
#' @return contours, oriented correctly
#'
enforce_output_contour_orientation<-function(contours, fault_trace){
    new_contours = contours
    fault_trace_coords = coordinates(fault_trace)[[1]]

    # Check that first unit vector of fault_trace_contours is oriented the same
    # as unit vector of new_contours
    n_ftc = fault_trace_coords[[1]][1,] - fault_trace_coords[[1]][2,]

    for(i in 1:length(contours)){
        contour_i_coords = coordinates(contours)[[i]][[1]]
        l = length(contour_i_coords[,1])
        n_c = contour_i_coords[1,] - contour_i_coords[2,]
        # If dot product is positive, they are oriented ok
        if(sum(n_c*n_ftc)<0.){
            new_contours@lines[[1]]@Lines[[1]]@coords = 
                contour_i_coords[l:1,]
        }
    }

    return(new_contours)
}

###############################################################################

option_list = list(

    make_option(
        c("-i", "--input_shapefile"), 
        type="character", 
        default=NA,
        help=paste0("Name of input shapefile defining the fault trace ",
                    "[default %default]")
    ),
    
    make_option(
        c("-d", "--fault_trace_depth"), 
        type="double", 
        default=0,
        help="Depth of fault trace in km [default %default]"
    ),


    make_option(
        c("-m", "--contour_min_depth"),
        type="double",
        default=5,
        help=paste0("Minimum output contour depth in km ( must be > ",
                    "fault_trace_depth) [default %default]")
    ),

    make_option(
        c("-M", "--contour_max_depth"),
        type="double", 
        default=60,
        help="Maximum output contour depth in km [default %default]"
    ),

    make_option(
        c("-c", "--contour_interval"), 
        type="double", 
        default=5,
        help=paste0("Output fault contour interval in km ",
                    "(approximate only, contours are forced to be evenly ",
                    "spaced between contour_min_depth and contour_max_depth) ",
                    "[default %default]")
    ),
    
    make_option(
        c("-o", "--output_shapefile"), 
        type="character", 
        default='trace_contours',
        help="Directory (and name) for output shapefile [default %default]"
    ),

    make_option(
        c("-v", "--verbose"),
        action="store_true",
        default=FALSE,
        help="Print extra output [default %default]"
    ),

    make_option(
        c("-p", "--plot_delay"), 
        type="double", 
        default=5.0,
        help="Seconds to display plot (if --verbose ) [default %default]"
    )

)

# Parse options
opt = parse_args(OptionParser(option_list = option_list))

fault_trace_shapefile = opt$input_shapefile 
fault_trace_depth = opt$fault_trace_depth
contour_interval = opt$contour_interval
contour_maxdepth = opt$contour_max_depth
contour_mindepth = opt$contour_min_depth
verbose = opt$verbose
plot_delay = opt$plot_delay
output_shapefile = opt$output_shapefile

# Error checks
if(is.na(fault_trace_shapefile)){
    stop('Must provide an input shapefile')
}

if(!file.exists(fault_trace_shapefile)){
    stop(paste0("Cannot find input_shapefile ", fault_trace_shapefile, 
                "\n Check that this file exists")
        )
}

if (verbose) print("Reading shapefile ...")

fault_trace = readOGR(dsn=fault_trace_shapefile, 
    layer=gsub('.shp', '', basename(fault_trace_shapefile)),
    verbose=FALSE)

fault_trace_extent = extent(fault_trace)

in_0_360 = ( (fault_trace_extent@xmin >= 0) & 
             (fault_trace_extent@xmax <= 360) )

if ((fault_trace_extent@xmin < -180) | (fault_trace_extent@xmax > 360)){
    stop("Fault trace longitude must either be within [-180,180], or [0, 360]")
}

ftc = coordinates(fault_trace)

for (i in 1:length(ftc)){

    if (length(ftc[[i]])>2){
    stop(paste0("Fault trace segments must only be defined by 2 points",
                "Segment ", i, " appears to have ", length(ftc[[i]])))
    }

}

# Call main routine
output_contours_info = extend_trace_to_depth_contours(
    fault_trace, 
    maxdepth = contour_maxdepth + contour_interval, 
    mindepth = fault_trace_depth, 
    contour_maxdepth = contour_maxdepth,
    contour_mindepth = contour_mindepth,
    in_0_360 = in_0_360,
    contour_interval = contour_interval, 
    verbose = verbose,
    plot_delay = plot_delay)


if(verbose) print('Ensuring correct contour orientation')

output_contours = enforce_output_contour_orientation(output_contours_info[[3]],
     fault_trace)

if(verbose) print('Writing output shapefile ...')

writeOGR(output_contours, dsn=output_shapefile, 
         layer=output_shapefile, driver='ESRI Shapefile', overwrite=TRUE)

