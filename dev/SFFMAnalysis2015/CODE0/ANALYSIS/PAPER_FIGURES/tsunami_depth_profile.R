##################################################################
# Parameters from python

depthFun<-function(x){
    # Here we code up the same depth function as in python,
    # Except viewed from the perspective of 'base of trench has x = 0'
    # x going towards shore = positive
    # This is the orientation of the deformation data
    #
    #
    beach_slope=-1.0e-02
    beach_end=5000.
    shelf_start=0. # Shelf starts this far beyond the beach end
    shelf_start_depth=0. # Shelf starts at this depth 
    shelf_finish=200000. # Shelf ends this far beyond the beach end
    shelf_bottom=2000. # Shelf drops this much
    slope_finish=300000. # Slope ends this far beyond the beach
    slope_bottom=6000. # Slope drops this much
        
    trench_slope=(slope_bottom-shelf_bottom)/(slope_finish-shelf_finish)
    
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

