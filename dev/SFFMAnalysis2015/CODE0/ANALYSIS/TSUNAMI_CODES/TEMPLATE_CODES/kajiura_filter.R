######################################################################
#
# Here we have functions for nearest-neighbour and triangular interpolation
#
# FIXME: Can later extract these into a single function
#

nearest_neighbour_interpolation<-function(xy, vals, newPts){
    ## Function for neaest neighbour interpolation
    ## Allows me to do this cleanly
    ##
    ## xy = points (typically nx2 matrix)
    ## vals = values at xy (can be a matix with 1 or more colums)
    ## newPts = points where we want interpolated values (typically nx2 matrix)
    library(SearchTrees)

    xyTree=createTree(xy)
    newInds=knnLookup(xyTree, newx=newPts[,1], newy=newPts[,2], k=1)

    if(is.null(dim(vals))){
        return(vals[newInds])
    }else{
        return(vals[newInds,])
    }
}

triangular_interpolation<-function(xy, vals, newPts, useNearestNeighbour=TRUE){
    # Here we do triangular interpolation with delaunay triangulation
    # Unfortunately the look-up is slow. 
    # So as an alternative, find the 3 nearest neighbours of each point, and
    # interpolate from that (with limiting)
   
    if(dim(xy)[1]<3 | is.null(dim(newPts))){
        stop('Need at least 3 input xy points, and newPts should be a matrix')
    }

    if(is.null(dim(vals))){
        if(length(vals)!=length(xy[,1])){
            stop('Length of xy[,1] and vals must be the same')
        }
    }else{
        if(length(vals[,1])!=length(xy[,1])){
            stop('Length of xy[,1] and vals[,1] must be the same')
        }
        
    }
 
    # Flag to say whether we use nearest-neighbours to define the triangulation
    #nnSort=(useNearestNeighbour & length(xy[,1])>3)
    nnSort=useNearestNeighbour

    if(!nnSort){
        # Use geometry package triangulation
        # This is slow for large problems [tsearch is presently slow], but the
        # best approach
        library(geometry) 
        triIndices=delaunayn(xy)
        # Use barycentric coordinates from tsearch
        triOn=tsearch(xy[,1], xy[,2], triIndices, newPts[,1], newPts[,2], bary=TRUE)
        if(is.null(dim(vals))){
            # Vals is a vector
            final=vals[triIndices[triOn$idx,1]]*triOn$p[,1]+
                  vals[triIndices[triOn$idx,2]]*triOn$p[,2]+
                  vals[triIndices[triOn$idx,3]]*triOn$p[,3]

        }else{
            # Vals is a matrix
            final=matrix(NA,ncol=ncol(vals), nrow=nrow(newPts))
            for(i in 1:ncol(final)){
                final[,i]=vals[triIndices[triOn$idx,1],i]*triOn$p[,1]+
                          vals[triIndices[triOn$idx,2],i]*triOn$p[,2]+
                          vals[triIndices[triOn$idx,3],i]*triOn$p[,3]
            }
        }

    }else{
        # Hack to try to speed-up tsearch on large problems.
        # Instead of using t-search, find the 3 nearest neighbours and make a triangle
        #browser()
        library(SearchTrees)
        #triXY=(xy[triIndices[,1],1:2]+ xy[triIndices[,2],1:2] + xy[triIndices[,3],1:2])/3
        triTree=createTree(xy)
        # Lookup the nearest index on the tree
        lookupInds=knnLookup(triTree, newx=newPts[,1],newy=newPts[,2],k=3)

        ## Interpolate. 
        ## Get vertices
        p1=matrix(xy[lookupInds[,1],], ncol=2)
        p2=matrix(xy[lookupInds[,2],],ncol=2)
        p3=matrix(xy[lookupInds[,3],], ncol=2)
       
        ## Get triangle gradient  
        dx31=p3[,1]-p1[,1]
        dy31=p3[,2]-p1[,2]
        dx21=p2[,1]-p1[,1]
        dy21=p2[,2]-p1[,2]

        dxN=newPts[,1]-p1[,1]
        dyN=newPts[,2]-p1[,2]
        #
        ## Compute triangle area
        area= dx21*dy31-dx31*dy21 #dy21*dx31 - dx21*dy31

        # Gradient coefficients
        a = (dy31*dxN-dx31*dyN)/area 
        b = (-dy21*dxN+dx21*dyN)/area 

        # Treat cases with degenerate triangles -- use nearest-neighbour instead
        EPS=1.0e-06
        a[abs(area)<EPS]=0.
        b[abs(area)<EPS]=0.


        ##b = (dx21*dz31 - dx31*dz21)/area
        if(is.null(dim(vals))){
            # Find max/min 'vals' on triangle
            valsMax=pmax(vals[lookupInds[,1]], vals[lookupInds[,2]], vals[lookupInds[,3]])
            valsmin=pmin(vals[lookupInds[,1]], vals[lookupInds[,2]], vals[lookupInds[,3]])
            
            dz31=vals[lookupInds[,3]] - vals[lookupInds[,1]]
            dz21=vals[lookupInds[,2]] - vals[lookupInds[,1]]

            final = vals[lookupInds[,1]] + a*dz21 + b*dz31

            # Limit
            M=final>valsMax
            m=final<valsmin
            limit=pmax(M,m)
            final = final*(1-limit) + vals[lookupInds[,1]]*limit
        }else{
            # Vals is higher dimensional

            # Find max/min 'vals' on triangle
            valsMax=pmax(vals[lookupInds[,1],], vals[lookupInds[,2],], vals[lookupInds[,3],])
            valsmin=pmin(vals[lookupInds[,1],], vals[lookupInds[,2],], vals[lookupInds[,3],])

            dz31=vals[lookupInds[,3],] - vals[lookupInds[,1],]
            dz21=vals[lookupInds[,2],] - vals[lookupInds[,1],]

            final = vals[lookupInds[,1],] + a*dz21 + b*dz31

            # Limit
            M=final>valsMax
            m=final<valsmin
            limit=pmax(M,m)
            # If outside min/max, use nearest neighbour only
            final = final*(1-limit) + vals[lookupInds[,1],]*limit

        #    return(final)
        }
    }
    return(final)
}

test_triangular_interpolation<-function(){
    # Make a single triangle in the plane z=x+y, and interpolate from it
    xy=matrix(c(0,0,0,1,1,1),ncol=2,byrow=T)
    vals=c(0, 1, 2) # z=x+y
    newPts=matrix(c(0.5, 0.5, 0.3, 0.3), ncol=2, byrow=T)

    out=triangular_interpolation(xy, vals, newPts)
    if(all(all.equal(out, c(1.0,0.6)))){
        print('PASS')
    }else{
        print(out)
        print('FAIL')
    }

    # Re-order triangle
    xy=xy[3:1,]
    vals=vals[3:1]
    out=triangular_interpolation(xy, vals, newPts)
    if(all(all.equal(out,c(1.0,0.6)))){
        print('PASS')
    }else{
        print(out)
        print('FAIL')
    }

    #another one, with formula z=0.5*x+0.2*y+7
    xy=matrix(c(-1, -1, 1, -0.5, 0.5, 1), ncol=2,byrow=2)
    vals=0.5*xy[,1]+0.2*xy[,2]+7
    newPts=matrix(c(0,0, 0.5, 0.3),ncol=2,byrow=T)
    expectedVals=0.5*newPts[,1]+0.2*newPts[,2]+7
    out=triangular_interpolation(xy,vals,newPts)
    if(all(all.equal(out,expectedVals))){
        print('PASS')
    }else{
        print(out)
        print('FAIL')
    }

    # A point outside the triangle 
    newPts=matrix(c(-1,0, -1, 1),ncol=2,byrow=T)
    out=triangular_interpolation(xy,vals,newPts)
    if(all(is.na(out))| formals(triangular_interpolation)$useNearestNeighbour){
        print('PASS')
    }else{
        print(out)
        print('FAIL')
        print('(  Failure is expected here if using approximate triangulation based on nearest neighbour methods)')
    }

    # A single point
    newPts=matrix(c(0,0),ncol=2)
    out=triangular_interpolation(xy,vals,newPts)
    if(out==7){
        print('PASS')
    }else{
        print(out)
        print('FAIL')
    }

    # Points on the triangle
    newPts=xy
    out=triangular_interpolation(xy,vals,newPts)
    if(all(out==vals)){
        print('PASS')
    }else{
        print('FAIL')
    }

    # Point on an edge
    newPts=matrix(0.5*(xy[1,]+xy[2,]),ncol=2)
    out=triangular_interpolation(xy,vals,newPts)
    if(all(out==0.5*(vals[1]+vals[2]))){
        print('PASS')
    }else{
        print(out)
        print('FAIL')
    }

}






#######################################################################
#
# Here we try to implement an initial water surface deformation filter
#


kajiura_g<-function(r, nt=1000, verbose=FALSE, recursive_stop_factor=1.0e-05){
    # Compute kajiura's filter function
    if(length(r)>1) stop('r can only be length 1')
    # 

    # The filter has an infinite series representation
    n=seq(0,nt)
    g=1/pi*sum( (-1)^n*(2*n+1)/( (2*n+1)^2 + r^2)^(3/2))

    n2=seq(nt+1,2*nt)
    g2=1/pi*sum( (-1)^n2*(2*n2+1)/( (2*n2+1)^2 + r^2)^(3/2))
        
    if(verbose) print(c(g2,g))

    if(abs(g2)>abs(g)*recursive_stop_factor){
        # Ensure that we have enough terms
        if(verbose) print('recursing')
        g=kajiura_g(r,nt=2*nt,verbose=verbose)
    }
    return(g+g2) 
} 

kg_empirical<-function(xMax=9, n=81){
    # Empirical approximation to kajiura_g
    #
    # Uses splines to log-transformed data
    #
    # It's good
    #
    kg = Vectorize(kajiura_g)
    x = seq(0,xMax,len=n)
    kg_x = kg(x)
    f_kx = splinefun(x, log(kg_x))
    kgE<-function(x) exp(f_kx(x))

    return(kgE)
}

test_kajiura_g<-function(){ 
    # Test that kajiura_g is correctly coded 
    testr=4 
    # (1/2pi)* (The integral of this) is another representation of the G function 
    f<-function(m) m*besselJ(m*testr,nu=0)/cosh(m) 
    fInt=integrate(f, 0, Inf) 
    fInt=fInt$value/(2*pi) 
    ref_val=kajiura_g(testr) 
    out=abs(fInt-ref_val) 
    print('Compare theoretical methods 1 and 2')
    print('Difference is: ')
    print(out)
    print('Differences < 1.0e-03 are ok')
    if(out<1.0e-03){
        print('PASS')
    }else{
        print('FAIL')
    }

    print('Compare theoretical and empirical')
    kg_E=kg_empirical()
    # Add small perturbation to x point so that we don't evaluate at spline knot
    e_ref=kg_E(testr+0.001)
    ref_val=kajiura_g(testr+0.001)
    out=abs(ref_val-e_ref) 
    print('Difference is: ')
    print(out)
    print('Differences < 1.0e-06 are ok')
    if(out<1.0e-06){
        print('PASS')
    }else{
        print('FAIL')
    }

    x=seq(0,8,len=50)
    kgV=Vectorize(kajiura_g)
    theory=kgV(x)
    splineApprox=kg_E(x)
    par(mfrow=c(2,2))
    plot(x,theory,t='o',log='y',main='Kajuira function (log y-scale)')
    points(x,splineApprox,t='l',col=2)
    legend('topright',c('Theory', 'Splines'), lty=c(1,1),col=c(1,2)) 

    plot(theory-splineApprox, main='ABS Difference',t='l')
    plot((theory-splineApprox)/theory, main='REL Difference',t='l')
} 

#######################################################
#
#
# MAIN FUNCTION BELOW
#
#
#######################################################

tsunami_init_filter<-function(xyDef, 
                              depth,
                              grid_dx=max(depth)/2,
                              grid_dy=max(depth)/2, 
                              edge_buffer_value=0, 
                              edge_effect_correction_scale=1.5, 
                              kajiuraGmax=9,
                              interpolator='linear',
                              volume_change_error_threshold=0.02,
                              verbose=TRUE){ #){ 
    #
    # Implement the filter similar to that of Glimsdal et al (2013), based on kajiura (1963)
    #
    # newDeformation(x,y) = int int (oldDeformation(x',y')*G( sqrt( (x'-x)^2+(y'-y)^2) / depth(x, y) ) dx' dy'
    # 
    # This is a 2D generalisation of the cosh filter, justified for a 'temporally short'
    # earthquake with ocean governed by linear shallow water equations in constant depth water
    #
    # Essentially:
    #    xyDef[,3] <-- convolution of [ (xyDef[,3]) and G(r/depth) ] where G is
    #    a filter function (kajiura's G) adjusted to have integral 1 over the
    #    filter window
    #
    # We attempt to reduce edge effects with linearly weighting old and filtered values at edges
    #  since we cannot EFFICIENTLY deal with edge effects in a better way. Therefore
    # it is best to have unimportant features around the edge of the input points
    #
    # We actually allow xyDef to be unstructured, and start by gridding the results
    # on a grid with spacing approximately grid_dx,grid_dy, 
    #  , using nearest neighbour interpolation. Generally, setting these values
    #  to a fraction of the refernece depth should be ok. Something a bit
    #  smaller than the input resolution would be good.
    # The grid spacing is not exactly grid_dx,grid_dy, because it is forced to be 
    #   an integer divisor of reference_depth
    #
    # NOTE that because we re-grid the deformation/depth data, we assume that both are continuous
    # 
    # INPUTS:
    # xyDef  = 3 column matrix with x,y, Deformation. 
    #          Assumed to lie on a grid
    #
    # depth = vector with the the depth at each x,y point in xyDef
    #
    #
    # grid_dx, grid_dy = Point spacing in the 'internal' gridded data used for smoothing
    #
    #
    # edge_buffer_value -- Outside the domain edges we assume this is the value
    #              of xyDef, when the filter is applied. Without some
    #              correction, this would make edges tend towards
    #              edge_buffer_value
    #
    # edge_effect_correction_scale -- To reduce edge effects, after filtering, we compute 
    #           the distance of every point to the edge of the domain 'd', and 
    #           then return the solution 
    #           WT = max( 1-d/(reference_depth*edge_effect_correction_scale), 0)**0.5
    #           OUTPUT = WT*(ORIGINAL OUTPUT) + (1-WT)*FILTERED OUTPUT
    #
    # interpolator -- 'linear' or 'nearest'. Linear is better, but slow for
    #                  large point clouds. Don't use nearest unless grid_x,
    #                  grid_y are 'small enough'
    #
    # kajiuraGmax -- When empirically approximating kajiuraG, we fit it from
    #                  x=[0, kajiuraGmax]. Values above this are evaluated to zero
    #
    # volume_change_error_threshold -- If the difference in the volume before
    #            and after filtering, relative to the original volume, is less than this,
    #            then throw an error
    #
    # verbose == print lots of stuff
    # 
    # OUTPUTS: replacement version of xyDef, with smoothing applied


    #browser()

    reference_depth=max(depth)

    # Get fast approximation to G function
    kgE=kg_empirical(xMax=kajiuraGmax)

    # dx/dy for gridded data + filter
    dx=reference_depth/ceiling(reference_depth/grid_dx)
    dy=reference_depth/ceiling(reference_depth/grid_dy)

    m0=min(xyDef[,1])
    m1=max(xyDef[,1])
    newX=seq(m0,m1,len=round((m1-m0)/dx)+1)
    m0=min(xyDef[,2])
    m1=max(xyDef[,2])
    newY=seq(m0,m1,len=round((m1-m0)/dy)+1)
    lny=length(newY)
    lnx=length(newX)

    # Compute nearest-neighbour interpolation
    newPts=as.matrix(expand.grid(newX,newY))

    # Create search tree and get values on new grid
    if(verbose) print('Unstructured interpolation number 1...')
    if(interpolator=='nearest'){
        # Nearest neighbour interpolation
        interp1=nearest_neighbour_interpolation(xyDef[,1:2], cbind(xyDef[,3], depth), newPts)
        newVals=matrix(interp1[,1], ncol=lnx,byrow=T)
        newDepth=matrix(interp1[,2], ncol=lnx,byrow=T)
    }else if(interpolator=='linear'){
        # Delaunay triangulation interpolation
        interp1=triangular_interpolation(xyDef[,1:2], cbind(xyDef[,3], depth), newPts)
        # Remove NA values
        interp1[is.na(interp1[,1]), 1] = edge_buffer_value
        interp1[is.na(interp1[,2]), 2] = 1.e-12 # Set depth to a small number

        newVals=matrix(interp1[,1], ncol=lnx,byrow=T)
        newDepth=matrix(interp1[,2], ncol=lnx,byrow=T)
    }else{
        stop('interpolator not recognized')
    }
   
    # Compute kajiura_g filter function, on a matrix varying from +- 6 reference depths
    # Glimsdal et al highlight that this only needs to be 5 reference depths long
    # I agree from my own checks, but use 6 for safety [why? FIXME: It's substanially slower]
    filter_refDepth_range = 6
    fR = filter_refDepth_range*reference_depth
    filterXs = seq(-fR,fR, by=dx)
    filterYs = seq(-fR,fR, by=dy)
    lfx = length(filterXs)
    lfy = length(filterYs)

    if((length(filterXs)%%2 != 1)|
       (length(filterYs)%%2 != 1)){
        stop('ERROR: Filter size is wrong')
    }

    # Compute 'radius' term on filter
    filterXY = expand.grid(filterXs,filterYs)
    filterXYr = matrix( (filterXY[,1]^2+filterXY[,2]^2)**0.5,
                      ncol=lfx, byrow=TRUE)


    if(verbose) print('Applying filter ...') 
    # To apply the filter, we 'pad' the data with zeros (which means the edges
    # will taper down). 
    # The matrix has length(filterYs) rows of zeros at the top an bottom, and 
    #    length(filterXs) columns of zeros at the left and right
    filVals = matrix(edge_buffer_value, ncol=lnx+2*lfx, nrow=lny+2*lfy)
    filVals[lfy+1:lny, lfx+1:lnx] = newVals


    old_newVals = newVals

    # Now set newVals to zero -- it will hold the filtered results
    newVals = 0.*newVals
    GtermsSum = 0.*newVals

    # The filter is of length lfx,lfy
    # Generally lfx<<lnx, lfy<<lny
    # So for efficiency, here we loop over every element of the filter when
    # computing the weighted average
    for(i in 1:lfx){
        if(verbose) print(paste(i, ' of', lfx ))
        for(j in 1:lfy){
            # Compute r/depth for the j,i cell of the filter,  avoid division by zero
            r_on_d = pmin( filterXYr[j,i]/pmax(newDepth,1.0e-20), kajiuraGmax)

            # Put into a matrix which aligns with newVals
            G_j_i = matrix(kgE(r_on_d), ncol=lnx)

            # Numerator of the weighted average
            newVals = newVals + filVals[(j+(lfy-1)/2)+1:lny,(i+(lfx-1)/2)+1:lnx]*G_j_i

            # Denominator of the weighted average. 
            GtermsSum = GtermsSum + G_j_i
        }
    }

    # Compute final weighted average 
    newVals = newVals/(GtermsSum)

    if(edge_effect_correction_scale>0.){
        if(verbose) print('Reducing edge effects ...')
        # Reduce edge-effects with a weighted average of the old values there
        # Compute the distance from an edge in the x/y directions
        xEdge = matrix( pmin((1:lnx)-0.5, (lnx:1)-0.5)*dx,byrow=T,ncol=lnx,nrow=lny)
        yEdge = matrix(pmin((1:lny)-0.5, (lny:1)-0.5)*dy,byrow=F,ncol=lnx,nrow=lny)
        # Convert to a weight
        xEdge = pmax(1-xEdge/(fR*edge_effect_correction_scale),0.)**0.5
        yEdge = pmax(1-yEdge/(fR*edge_effect_correction_scale),0.)**0.5
        edgeF = pmax(xEdge,yEdge)
        # Take weighted average of the original values and the new ones
        newVals = edgeF*old_newVals+(1-edgeF)*newVals
    }
   

    # Back-compute x,y,Def with nearest neighbour lookup
    if(verbose) print('Unstructured interpolation number 2...')
    new_xyDef = xyDef
    #newXY_Tree=createTree(newPts)
    #new_xyDef_inds=knnLookup(newXY_Tree,newx=xyDef[,1],newy=xyDef[,2],k=1)
    if(interpolator=='nearest'){
        interp2 = nearest_neighbour_interpolation(newPts, c(t(newVals)), xyDef[,1:2])

        new_xyDef[,3] = interp2 #c(t(newVals))[new_xyDef_inds]
    }else if(interpolator == 'linear'){
        interp2 = triangular_interpolation(newPts, c(t(newVals)), xyDef[,1:2])
        interp2[is.na(interp2)] = edge_buffer_value
        new_xyDef[,3] = interp2 #c(t(newVals))[new_xyDef_inds]
    }

    newValsPosSum = sum(newVals*(newVals>0))
    newValsNegSum = sum(newVals*(newVals<0))
    old_newValsPosSum = sum(old_newVals*(old_newVals>0))
    old_newValsNegSum = sum(old_newVals*(old_newVals<0))

    oldNewvalsSum = sum(old_newVals)
    newvalsSum = sum(newVals)

    r1 = (newValsPosSum-old_newValsPosSum)/old_newValsPosSum
    r2 = (newValsNegSum-old_newValsNegSum)/old_newValsNegSum
    r3 = (sum(newVals)-sum(old_newVals))/sum(old_newVals)

    if(verbose) print(paste('Original re-gridded volume: ', oldNewvalsSum))
    if(verbose) print(paste('New gridded volume: ', newvalsSum))
    if(verbose) print(paste('Total volume relative error: ',r3))
    if(verbose) print(' ')
    if(verbose) print(paste('Original positive re-gridded volume: ', old_newValsPosSum))
    if(verbose) print(paste('New positive gridded volume: ', newValsPosSum))
    if(verbose) print(paste('Positive volume relative error: ', r1))
    if(verbose) print(' ')
    if(verbose) print(paste('Original negative re-gridded volume: ', old_newValsNegSum))
    if(verbose) print(paste('New negative gridded volume: ', newValsNegSum))
    if(verbose) print(paste('Negative volume relative error: ', r2))
    if(verbose) print(' ')

    # Check for gross volume conservation errors
    # In general, we cannot exclude this to some extent
    if(is.finite(r3)){
        if(abs(r3)>volume_change_error_threshold){
            print(paste('r3 = ', r3))
            stop('Volume change error threshold exceeded')
        }
    }else{
        # Look for other measures of problems
        # Note though that smoothing can decrease both the positive
        # and negative volume
        if(is.finite(r1)){
            if(abs(r1)>volume_change_error_threshold){
                print(paste('r1 = ', r1))
                stop('(Positive vol) Volume change error threshold exceeded')

            }
        }
        if(is.finite(r2)){
            if(abs(r2)>volume_change_error_threshold){
                print(paste('r2 = ', r1))
                stop('(Negative vol) Volume change error threshold exceeded')

            }
        }

    }

    return(new_xyDef)
     
}

##

check_tsunami_init_filter<-function(myFile="DEFORMATION/139944854082882_NST_clip_SSD_FALSE_RCS_TRUE/S_1923009001_M7.95_Kanto-Japan_65/OceanInitial_1.xyz"){
    # Basic code to check that the filter function qualitatively does what it
    # should
    xyDef=matrix(scan(myFile,sep=","),ncol=3,byrow=T)

    m1=min(xyDef[,2])
    m2=max(xyDef[,2])

    # Make depth vary with y
    var1=(xyDef[,2]-m1)/(m2-m1)
    depth=var1*(-3200)+3000 -300*sin(var1*2*pi*5)
    depth=pmax(depth,0)

    newxyDef=tsunami_init_filter(xyDef, depth,verbose=T, grid_dx=1500, grid_dy=1500)

    m1=matrix(xyDef[,3],ncol=500)
    m2=matrix(newxyDef[,3],ncol=500)
  
    print('Range of initial data is') 
    print(range(m1))
    print('Range of new data is')
    print(range(m2))

    print('Range of difference is')
    print(range(m2-m1))

    plot(m1[250,],t='l')
    points(m2[250,],t='l',col=2)

    return(environment())
 
}

test_tsunami_init_filter_2<-function(){
    #
    # Smoothing of step function
    #
    #

    lScale=50000
    n=100
    pts=as.matrix(expand.grid(seq(0,1,len=n), seq(0,2,len=n)))*lScale
    pts_z=pts[,1]>lScale/2.
    xyDef=cbind(pts,pts_z)
    
    # Make the depth constant at the top / bottom, and otherwise
    # linearly varying
    depth=200+pmin( pmax(xyDef[,2]-0.2*lScale, 0), 1.6*lScale)/100

    new_xyDef=tsunami_init_filter(xyDef, depth, grid_dx=400,grid_dy=300)
    m1=matrix(new_xyDef[,3],ncol=n)

    new_xyDef=tsunami_init_filter(xyDef, depth, grid_dx=200,grid_dy=200)
    m2=matrix(new_xyDef[,3],ncol=n)
    
    depthMat=matrix(depth,ncol=n)
    #new_xyDef=tsunami_init_filter(xyDef, depth, grid_dx=100,grid_dy=100)
    #m3=matrix(new_xyDef[,3],ncol=100)

    par(mfrow=c(2,1))    
    plot(m1[,50],t='l')
    points(m2[,50],t='l',col=2)
    #points(m3[,50],t='l',col=3)    

    image(m2,main='Water deepens north so edge gets smoother (but note edge effects)')

    # Put a check to see if they are the 'same' away from edges
    testStat=max(abs(m1[,20:80]-m2[,20:80]))
    if(testStat>0.1){
        print(testStat)
        print('FAIL')
    }else{
        print('PASS')
    }

    return(environment())
}


test_tsunami_init_filter_1<-function(){
    # Check that the smoothing of a 'point' source    
    # is as desired, with an uneven point cloud
    #
    # NOTE: The volume of the 'point' source is not well defined.
    #       Since the input xyDef data is interpreted as unstructured points,
    #        the volume of a 'point' source depends on how you interpolate the data
    #       So, we account for this in the test
    #
    # However, the repeated re-gridding in the routine will still introduce some errors

    # DEFINE INITIAL xyDef + depth
    xrange=seq(-1,1,len=201)*50000
    yrange=seq(-1,1,len=401)*50000

    xyGrid=as.matrix(expand.grid(xrange,yrange))

    mm=which((xyGrid[,1]==0) & (xyGrid[,2]==0))
    if(length(mm)!=1) stop('Did not find origin')

    # Make deformation mostly zero, except for a 'spike' of 1 at 0,0
    depths=xyGrid[,1]*0+2000
    def=xyGrid[,1]*0
    def[mm]=1.0

    xyDef=cbind(as.matrix(xyGrid), def)

    # Interpolate / smooth
    # NOTE: We need to ensure grid_dx, grid_dy are small enough to 'see'
    # the input data correctly. 
    new_xyDef=tsunami_init_filter(xyDef, depths, grid_dx=500,grid_dy=250)

    # Unsmoothed
    m1=matrix(xyDef[,3],ncol=201,byrow=T)
    # Smoothed
    m2=matrix(new_xyDef[,3],ncol=201,byrow=T)
    #
    d1=matrix(depths, ncol=201,byrow=T)

    # Find the 'r' value to apply to kajiura G
    rMax=matrix(sqrt(xyGrid[,1]^2+xyGrid[,2]^2), ncol=201,byrow=T)/d1

    # Make kajiura function
    kgE=kg_empirical()

    # Limit the bounds of the kajiura function inputs
    rMax=pmin(rMax,8)

    k2=matrix(kgE(rMax)*(rMax<8),ncol=201)
    k2=k2/sum(k2)*sum(m2) # Normalise to sum, to get around the 'volume' interpretation issue for a point source

    testStat=max(abs(m2-k2))
    if(testStat<1.0e-05){
        print('PASS')
    }else{
        print(testStat)
        print('FAIL')
    }

    return(environment())
}

test_tsunami_init_filter_3<-function(){
    # Check that the smoothing of a 'point' source    
    # is as desired, with an uneven point cloud
    # USE VARIABLE DEPTHS
    #
    # NOTE: The volume of the 'point' source is not well defined.
    #       Since the input xyDef data is interpreted as unstructured points,
    #        the volume of a 'point' source depends on how you interpolate the data
    #       So, we account for this in the test, by ensuring sum(k2)=sum(m2)
    #
    # However, the repeated re-gridding in the routine will still introduce some errors

    # DEFINE INITIAL xyDef + depth
    xrange=seq(-1,1,len=201)*50000
    yrange=seq(-1,1,len=401)*50000

    xyGrid=as.matrix(expand.grid(xrange,yrange))

    mm=which((xyGrid[,1]==0) & (xyGrid[,2]==0))
    if(length(mm)!=1) stop('Did not find origin')

    # Make deformation mostly zero, except for a 'spike' of 1 at 0,0
    def=xyGrid[,1]*0
    def[mm]=1.0
   
    # Make 'varying' depths 
    # If they vary too rapidly, the volume error becomes too large
    #browser()
    depths=2000.+100.*(sin(xyGrid[,1]/3000)+cos(xyGrid[,2]/3000*0.8)) #runif(length(xyGrid[,1]), min=-1,max=1)*500.+2000
    xyDef=cbind(xyGrid, def)

    # Interpolate / smooth
    # NOTE: We need to ensure grid_dx, grid_dy are small enough to 'see'
    # the input data correctly. 
    new_xyDef=tsunami_init_filter(xyDef, depths, grid_dx=500,grid_dy=250)

    # Unsmoothed
    m1=matrix(xyDef[,3],ncol=201,byrow=T)
    # Smoothed
    m2=matrix(new_xyDef[,3],ncol=201,byrow=T)
    #
    #d1=matrix(depths, ncol=201,byrow=T)
    #
    # Find the 'r' value to apply to kajiura G
    #rMax=matrix(sqrt(xyGrid[,1]^2+xyGrid[,2]^2), ncol=201,byrow=T)/d1
    rMax=sqrt(xyGrid[,1]^2+xyGrid[,2]^2)/depths
    # Limit the bounds of the kajiura function inputs
    rMax=pmin(rMax,8)

    # Make kajiura function
    kgE=kg_empirical()


    k2=matrix(kgE(rMax)*(rMax<8), ncol=201,byrow=T)
    k2=k2/max(k2)*max(m2) # Normalise to sum, to get around the 'volume' interpretation issue for a point source
    
    # Add a graphic
    par(mfrow=c(2,2))
    plot(m2[200,],t='l', main='A transect')
    points(k2[200,],t='l',col=2)
    plot(m2[,100],t='l', main='Another transect')
    points(k2[,100],t='l',col=2)
    mMax=max(max(m2), max(k2))
    image(m2, zlim=c(0,mMax), asp=1)
    image(k2, zlim=c(0,mMax), asp=1)

    # Test = difference in surfaces / maximum.
    testStat=max(abs(m2-k2))/max(m2)
    if(testStat<1.0e-01){
        print('PASS')
    }else{
        print(testStat)
        print('FAIL')
    }


    return(environment())
}


testAll<-function(){
    pdf('testAll.pdf',width=12,height=8)
    print(' ')
    print('Testing Triangular interpolation')
    test_triangular_interpolation()
    print(' ')
    print('Testing G function computation')
    test_kajiura_g()
    print(' ')
    print('Testing filter_3')
    T3=test_tsunami_init_filter_3()
    print(' ')
    print('Testing filter_2')
    T2=test_tsunami_init_filter_2()
    print(' ')
    print('Testing filter_1')
    T1=test_tsunami_init_filter_1()
    print(' ')
    print('Running a realistic case, if it fails, that is an error...')
    C1=check_tsunami_init_filter()
    dev.off()
    return(environment())
}


