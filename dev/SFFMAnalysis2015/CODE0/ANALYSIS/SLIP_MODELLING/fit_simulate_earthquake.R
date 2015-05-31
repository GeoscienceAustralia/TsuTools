###############################################################################
#'
#' Code to fit/simulate spectral models of earthquake slip surfaces
#' , and wrappers to the Okada function
#'
#' AUTHORS
#' Gareth Davies, Geoscience Australia 2013-2015, gareth.davies.ga.code@gmail.com
#'
#' @@ other contributors add names above here @@
#'

# FIXME: Try passing all packages to cluster with clusterEvalQ
suppressPackageStartupMessages(library(rgdal))
suppressPackageStartupMessages(library(raster))
suppressPackageStartupMessages(library(stabledist))

if (!exists('config_pars')){
    stop('config_pars must exist for fit_simulate_earthquake.R to be used')
}

###############################################################################
#'
#' Utility function to get the random seed -- and make it first if required
#'
get_random_seed<-function(){

    if(!exists('.Random.seed', where=.GlobalEnv)){
        # Create the random seed by calling a random number generator
        x = rnorm(1)
    }

    return(.Random.seed)
}


#' Utility function to change the random seed (by asking for a random number)
#'
#' This was useful with the mc-apply and co family of functions for parallel
#' computing (replaced with clusterMap etc)
#' Our code mainly uses random numbers inside parallel
#' loops -- so the master node didn't actually have its seed updated. 
#' However, the seeds on the worker nodes are set (uniqely) based on the master node seed.
#' As a result, the same random numbers are often reused in different calls to
#' mclapply. (Although within a single call to mclapply and similar, the random
#' generator uses different seeds on different processors -- so we still have
#' 'random' numbers within each call)
#'
#' For example, if I call:
#' mclapply(as.list(12), f<-function(x) mean(.Random.seed), mc.cores=12)
#' Then I get a list of length 12, and all the results are different. 
#' But if I call it again, I get the same result (since the .Random.seed on the master node is not changed)
#' On the other hand, this doesn't happen using 'clusterMap' et al.
#'
change_random_seed<-function(){

    x = rnorm(1)

    return(.Random.seed)
}


###############################################################################
#'
#' Compute wavenumbers in numerical space.
#' i.e. 0, 1/N, 2/N, ...
#'
#' Can be adjusted to physical wavenumbers by division by cellsize dx for
#' kcx (or dy for kcy) 
#
#' @param tg_mat = matrix or RasterLayer
#' @return list of wavenumbers in numerical space
#'
get_wavenumbers<-function(tg_mat){

    #B1 = nrow(tg_mat)
    #A1 = ncol(tg_mat)
    tmp = dim(tg_mat)
    B1 = tmp[1]
    A1 = tmp[2]
    X2 = pmin(0:(B1-1), B1-(0:(B1-1)))/B1
    X1 = pmin(0:(A1-1), A1-(0:(A1-1)))/A1
    kx = matrix(X1,ncol=A1,nrow=B1, byrow=T)
    ky = matrix(X2, ncol=A1,nrow=B1)

    return(list(kx,ky))

}

###############################################################################
#'
#' Make a random slip surface from a template raster, given some regression
#' parameters (specified in NUMERICAL SPACE independent of the pixel size) and 
#' SFFM definitions in config_pars
#'
#' If kx = 0,1/N,2/N, ... are the numerical wavenumbers,
#' they are equal to the ('physical' wavenumbers) x (delX)
#' In that case, if reg_par[1] and reg_par[2] are the numerical corner 
#' wavenumbers, then:
#' reg_par[1] = kcx*delX 
#' [where kcx is the physical 'corner wavenumber' and delX is the grid X size.]
#' And similarly for reg_par[2] with kxy/delY instead of kcx/delX
#'

#' @param reg_par is a vector of c(kcxN, kcyN), in NUMERICAL SPACE (or
#'        a 3rd parameter is accepted in some cases, see code)
#' @param tg_mat is a 'template' raster, or matrix
#' @param local_hurst value for the hurst exponent
#' @return Output is the same class as tg_mat
#'
simulate_sffm<-function(reg_par, tg_mat, local_hurst = config_pars$hurst){ 

    # Record random seed for reproducibility
    initial_seed = get_random_seed()
    
    # Make wavenumber matrices
    tmp = get_wavenumbers(tg_mat)
    kx = tmp[[1]]
    ky = tmp[[2]]


    # Compute modelled fourier spectrum
    model_fourierspec = config_pars$model_fourierfun(kx, ky, reg_par, 
        hurst = local_hurst)

    # Note: if kx = 0,1/N,2/N, ... are the numerical wavenumbers,
    # they are equal to the ('physical' wavenumbers) x (delX)
    # In that case, 
    # reg_par[1] = kx_Numerical = kcx*delX where kcx is the physical 'corner
    # wavenumber'

    switch(config_pars$noise_distribution,
        'gaussian' = {
            # Gaussian with mean = 0, sd = 1
            # Ultimately the sd doesn't matter for mu = 0
            # since changing the sd is equivalent to rescaling the random
            # sample, which doesn't effect the Arg(fft()) operation
            #r_noise = matrix(rnorm(length(tg_mat)), ncol = ncol(tg_mat), 
            #    nrow = nrow(tg_mat))
            r_noise = rnorm(length(tg_mat))
            dim(r_noise) = dim(tg_mat)[1:2]

            r_noise = r_noise*sign(sum(r_noise))
        },
        'cauchy' = {
            ## Parameters from Geist (2013)
            #cauchy_location = -1
            #cauchy_scale = 9.7

            # Parameters from median of ML fits to the data
            # If the location is 0., the scale is irrelevant since
            # the phase does not change under positive rescaling 
            cauchy_location = 0.
            cauchy_scale = 8.1

            #cauchy_random = rcauchy(length(tg_mat), location = cauchy_location,
            #    scale = cauchy_scale)
            #
            #r_noise = matrix( cauchy_random, ncol = ncol(tg_mat),
            #    nrow = nrow(tg_mat))
            r_noise = rcauchy(length(tg_mat), location = cauchy_location,
                scale = cauchy_scale)
            dim(r_noise) = dim(tg_mat)[1:2]

            r_noise = r_noise*sign(sum(r_noise))
        },
        'laplace' = {
            # Random laplace variables (can rescale later)
            #laplace_random = rexp(length(tg_mat))*
            #    sample(c(-1,1), size = length(tg_mat), replace=TRUE)
            #
            #r_noise = matrix(laplace_random, ncol = ncol(tg_mat),
            #    nrow = nrow(tg_mat))
            r_noise = rexp(length(tg_mat))*
                sample(c(-1,1), size = length(tg_mat), replace=TRUE)
            dim(r_noise) = dim(tg_mat)[1:2]

            r_noise = r_noise*sign(sum(r_noise))

        },
        'stable' = {
            # Random variables with stable distribution
            #
            # Here delta is a location parameter, and gamma is a scale parameter
            # With delta = 0, any gamma > 0 will give the same phase

            if(length(reg_par)>2){
                alpha_val = reg_par[3]
            }else{
                alpha_val = 1.81
            }
            #stable_random = rstable(length(tg_mat), alpha=reg_par[3], beta=0.0, 
            #    gamma=1.0, delta=0)
            #
            #r_noise = matrix(stable_random, ncol = ncol(tg_mat),
            #    nrow = nrow(tg_mat))
            r_noise = rstable(length(tg_mat), alpha=alpha_val, beta=0.0, 
                gamma=1.0, delta=0)
            dim(r_noise) = dim(tg_mat)[1:2]

            r_noise = r_noise*sign(sum(r_noise))

        },
        { stop('config_pars$noise_distribution not recognized')}
    )


    if(TRUE){
        # Use the random variables to generate the phase
        # This preserves the model spectral relation
        # (prior to clipping / filtering / etc)
        # This is the technique used by Gallovic, Mai, ...
        fake_phase = Arg(fft(r_noise))
        tmp = model_fourierspec*(cos(fake_phase)+1i*sin(fake_phase))
    }else{
        # Use the random variables and multiply in fourier space
        # This is the technique used by geist / Lavallee / ....
        # It can lead to a heavier skew in the slip (so truncation is normally
        # applied for stable distribution)
        tmp = model_fourierspec*fft(r_noise)
    }

    # Initial SFFM
    fake_data = Re(fft(tmp, inverse=T))/prod(dim(tg_mat))

    if(config_pars$spatial_slip_decay == 'inverse-box-cox'){
        # Box-cox type transformation must occur before negative value removal,
        # since the mean/sd matching steps (after Goda et al., 2014) can create
        # negative slip values
    
        # THEORY: 
        #
        # Recall y = exp(x) = lim(n-->Inf) (1+x/n)**n 
        #                   = lim(lambda-->0+) (1+lambda*x)**(1/lambda)
        #
        # This is not defined/sensible for x < -n, but we can do
        #    f = f(x) = (1+abs(x)/n)**(n*sign(x))
        # When x is negative, this evaluates to 1/f(abs(x)), as with the
        # exponential function. Here again n = 1/lambda. I've plotted this
        # function for various finite n, and it seems well behaved
        #
        # The box-cox transform (as reported often) is like the inverse of this
        # formula, with lambda = 1/n. 
        # x = (y**(lambda) -1) / lambda
        #
        # So as lambda-->0, the box-cox transform approximates the logarithmic
        # transform

        data_mean = mean(as.matrix(tg_mat))
        data_sd = sd(as.matrix(tg_mat))
        fake_data_coef_var = mean(fake_data)/sd(fake_data)

        # Set mean / sd from the data
        #fake_data = (fake_data-mean(fake_data))/
        #    sd(fake_data)*data_sd + data_mean
        #
        # If we set the sd from the data, then kcx/kcy seem to become
        # unidentifiable for some events -- and very very very small values
        # can be fit. This is because with very very small kcx/kcy, the
        # fake_data is nearly constant -- however, the rescaling by mean/sd 
        # undoes this.
        #
        # This alternative does not suffer that problem
        fake_data = fake_data - mean(fake_data)

        # Alternative: Exponential-like transformation, (inverse
        # of box-cox but applies to + and - numbers)
        lambda = config_pars$box_cox_lambda 
        fake_data = (lambda*abs(fake_data) + 1)**
            (sign(fake_data)/lambda)

        # Re-fix mean/sd 
        #fake_data = (fake_data-mean(fake_data))/
        #    sd(fake_data)*data_sd + data_mean
        
        # We don't in general know or model the sd -- so here we
        # rescale so the mean = data_mean, and coefficient of variation = 
        # initial_simulated_field coefficient of variation
        fake_data = (fake_data - mean(fake_data))/sd(fake_data)*
            (data_mean*fake_data_coef_var) + data_mean
    
        # There can still be negative values
    }


    # Deal with negative values
    switch(config_pars$negative_slip_removal_method,
        'abs' = {
            fake_data_clip = abs(fake_data)
        }, 
        'clip' = {
            # Clip negative values
            fake_data_clip = pmax(fake_data, 0)
        },
        'none' = {
            fake_data_clip = fake_data
        },
        stop('ERROR: config_pars$negative_slip_removal_method not recognized')
    )
  
    # Recentre the slip so that location of maxima in fake_data_clip is the
    # same as that in tg_mat.
    # This must happen after negative value removal, in case e.g. abs() leads
    # to a new location for the maxima
    if(config_pars$RECENTRE_SLIP){ 
        fake_data_clip = recentre_slip(fake_data_clip, tg_mat)
    }
    
    # Apply exponential slip decay. This must happen after recentering the slip
    if(config_pars$spatial_slip_decay %in% c('exponential', 'gaussian')){
        # Spatial filtering
        # Find distance from maxima
        maxInds = which(fake_data_clip == max(fake_data_clip), arr.ind=T)
        # Find distance from random point
        tmp = dim(fake_data_clip)[1:2]
        nc_fdc = tmp[2]
        nr_fdc = tmp[1]
        xMat = matrix(1:nc_fdc, ncol=nc_fdc, nrow=nr_fdc, byrow=TRUE)    
        yMat = matrix(1:nr_fdc, ncol=nc_fdc, nrow=nr_fdc, byrow=FALSE)    
        #
        # reg_par[1] = kcx*delx, so to get distance*kcx
        # [=distance/length_scale], we multiply the cell-count-distance by
        # reg_par[1].
        distMat2 = ( ((xMat-maxInds[2])*reg_par[1])**2 + 
                    ((yMat-maxInds[1])*reg_par[2])**2 )

        if(config_pars$spatial_slip_decay == 'exponential'){
            fake_data_clip = fake_data_clip*exp(-(distMat2**0.5))
        }else if(config_pars$spatial_slip_decay == 'gaussian'){
            fake_data_clip = fake_data_clip*exp(-(distMat2))
        }


    }else{
        # Just make sure config_pars has a valid value
        if(!(config_pars$spatial_slip_decay %in% c('none', 'inverse-box-cox'))){
            stop('config_pars$spatial_slip_decay not recognized')
        }
    }

    # Ensure final mean = data mean 
    if(class(tg_mat)=='RasterLayer'){
        fake_data_clip = fake_data_clip/sum(fake_data_clip)*sum(as.matrix(tg_mat))
    }else{
        fake_data_clip = fake_data_clip/sum(fake_data_clip)*sum(tg_mat)
    }

    if(class(tg_mat)=='RasterLayer' ){
        final_rast = raster(tg_mat)
        final_rast = setValues(final_rast, fake_data_clip)
    }else{
        final_rast = fake_data_clip
    }

    # Record the random number
    attr(final_rast, 'initial_seed') = initial_seed

    return(final_rast)        
}

#######################################################################
#'
#' Move asperities closer to the centre of the rupture,
#' or to the location of asperities in another slip raster (tg)
#'
#' If tg is not NULL, then move the max of m1 to the same location as
#' the max of tg, by re-ordering rows and columns
#' Otherwise, compute the max along each row -- then re-order the rows so that
#' the row with the smallest max is at the top of the rupture
#' (unless is is already at the top or bottom).
#' Then do the same with the columns, trying to put the smallest max
#' col on the left (unless it is already at the left or right)
#'
#' @param m1 rasterLayer or matrix containing a slip surface to be adjusted
#' @param tg rasterLayer or matrix containing a reference slip surface
#' @return recentred version of m1    
#'
recentre_slip<-function(m1, tg=NULL){

    if(!is.matrix(m1)){
        m1_mat=as.matrix(m1)
    }else{
        m1_mat = m1
    }

    if(is.null(tg)){
        # Find row with lowest max slip -- we want this on the bottom or top
        # edge
        row_to_bottom=which.min(apply(m1_mat,1,max))
        # Find col with lowest max slip -- we want this on the left or right
        # edge
        col_to_left=which.min(apply(m1_mat,2,max))
        nr=dim(m1_mat)[1]
        nc=dim(m1_mat)[2]
        if(!(row_to_bottom%in%c(1,nr))){
            #row_to_bottom is not on the bottom
            m1_mat=m1_mat[c( (row_to_bottom):nr, 1:(row_to_bottom-1)), ]
        }
        if(!(col_to_left%in%c(1,nc))){
            m1_mat=m1_mat[, c( (col_to_left):nc, 1:(col_to_left-1))]
        }

    }else{
        # Move the max to the location of the max of tg
        if(!is.matrix(tg)){
            tg_mat = as.matrix(tg)
        }else{
            tg_mat = tg
        }
        #new_row_max = which.max(apply(tg_mat,1,max)) 
        #new_col_max = which.max(apply(tg_mat,2,max)) 
        tmp = which(tg_mat == max(tg_mat), arr.ind=TRUE)
        new_row_max = tmp[1]
        new_col_max = tmp[2]

        #old_row_max = which.max(apply(m1_mat,1,max)) 
        #old_col_max = which.max(apply(m1_mat,2,max)) 
        tmp = which(m1_mat == max(m1_mat), arr.ind=TRUE)
        old_row_max = tmp[1]
        old_col_max = tmp[2]

        nr = dim(m1_mat)[1]
        nc = dim(m1_mat)[2]
        newRows = (((1:nr) - (new_row_max - old_row_max))%%nr)
        newRows[which(newRows==0)] = nr
        m1_mat = m1_mat[newRows,]
        newCols = (((1:nc) - (new_col_max-old_col_max))%%nc)
        newCols[which(newCols==0)] = nc 
        m1_mat = m1_mat[,newCols]
    }
    if(class(m1) == 'RasterLayer'){
        output = raster(m1)
        output = setValues(output,m1_mat)
    }else{
        # Work directly with matrices too
        output = m1_mat
    }
    return(output)
}

###############################################################################
#'
#' Function to get the along-strike / down-dip location of the peak slip
#' @param myrast raster of slip values
find_location_max<-function(myrast){
    # Find maxima
    maxinds = which(as.matrix(myrast)==max(as.matrix(myrast)),arr.ind=T)
    # Get x/y values in physical space -- use centre of pixel
    xval = (maxinds[2]-0.5)*res(myrast)[1]
    yval = (dim(myrast)[1] - maxinds[1]+0.5)*res(myrast)[2]
    return(c(xval,yval)) 
}


###############################################################################
#'
#' Plot slip vs distance from peak-slip location
#'
#'
slip_vs_distance<-function(myrast){

    slipMat = as.matrix(myrast)
    maxinds = which(slipMat == max(slipMat), arr.ind=TRUE)

    nx = dim(myrast)[2]
    ny = dim(myrast)[1]

    dx = res(myrast)[1]
    dy = res(myrast)[2]

    xMat = matrix(1:nx, ncol=nx, nrow=ny, byrow=TRUE)    
    yMat = matrix(1:ny, ncol=nx, nrow=ny, byrow=FALSE)    

    distMat = ( ((xMat-maxinds[2])*dx)**2 + 
                ((yMat-maxinds[1])*dy)**2 )**0.5

    return(list(distMat = distMat, slipMat = slipMat, xloc = xMat*dx, yloc = yMat*dy))
}

###############################################################################
#'
#' More useful stats
#'
Useful_eq_stats<-function(
    grid_files,
    asperity_scale=1.25, 
    asperity_cell_threshold=2, 
    recentre_FFM=FALSE){
    # Compute statistics about asperities / location of max slip / etc from
    # a vector of FFMs. 
    #
    # grid_files can either be a list of grid filenames, or a list of rasterLayers 

    Asp_area_store=rep(NA,length(grid_files)) # Area of 'asperity'
    total_80_area= rep(NA,length(grid_files)) # FFM area conttibuting 80% to the seismic moment
    Asperity_count=rep(NA,length(grid_files)) # Count asperities
    moment_scale=rep(NA,length(grid_files)) # Record stat related to seismic moment
    downDip_max_store=rep(NA,length(grid_files)) # distance of max slip DOWN FFM plane (not same as below ground)
    along_strike_max_store=rep(NA,length(grid_files)) # Along-strike location of max-slip 
    # Median slip in the highest-slip cells contributing 80% of the seismic moment
    median_80_slip=rep(NA,length(grid_files)) 
    # 95th percentile of slip in the highest-slip cells contributing 80% of the seismic moment
    quant_95_80_slip=rep(NA,length(grid_files)) 

    median_slip=rep(NA,length(grid_files))
    quant95_slip=rep(NA,length(grid_files))
    quant99_slip=rep(NA,length(grid_files))

    all_store=list() # Hold outputs
    for(i in 1:length(grid_files)){
        if(class(grid_files[i])=='character'){
            m1=raster(grid_files[[i]])
        }else{
            m1=grid_files[[i]]

        }

        if(recentre_FFM) m1=recentre_slip(m1)

        # Forcibly remove coordinate system -- can be wrong?
        proj4string(m1)=as.character(NA)


        # Find down-dip/along-strike location of max
        tmp=find_location_max(m1)
        downDip_max_store[i]=tmp[2]
        along_strike_max_store[i]=tmp[1]

        # Definition for asperities
        sorted_slip=sort(getValues(m1))
        cumsum_slip=cumsum(sorted_slip) # Sorted slip, cumulatively summed
        slip_80pc=sorted_slip[which(cumsum_slip>=0.2*max(cumsum_slip))]
        median_80_slip[i]=median(slip_80pc) #mean(slip_80pc)
        asperity_level=asperity_scale*median_80_slip[i] 

        # 95th percentile of slip on the 80%
        quant_95_80_slip[i] = quantile(slip_80pc, probs=0.95)

        median_slip[i]=quantile(sorted_slip,0.5)
        quant95_slip[i]=quantile(sorted_slip,0.95)
        quant99_slip[i]=quantile(sorted_slip,0.99)

        # Definition for moment scale = slip * area
        moment_scale[i]=sum(as.matrix(m1))*prod(res(m1))

        # Find regions>asperity_level
        m1_thresh=m1>asperity_level
        # Clump asperities (to isolate unique ones)
        m1_clump=as.matrix(clump(m1_thresh,gaps=FALSE))
        m1_clump[is.na(m1_clump)]=0 # Get rid of NA values

        # Compute area of asperity zones
        Asp_area_store[i]=sum(m1_clump>0,na.rm=T)*prod(res(m1)) 
        total_80_area[i]=length(slip_80pc)*prod(res(m1))
        m1_freq=table(m1_clump) # Statistics of distinct asperities

        # Include asperities if they have > asperity_cell_threshold cells
        tmp=m1_freq[2:length(m1_freq)]
        Asperity_count[i]=length(tmp[tmp>asperity_cell_threshold])
       
        # Convert asperity map back to raster 
        m1_clump2=raster(m1_thresh)
        values(m1_clump2)=m1_clump

        all_store[[i]]=list(m1=m1, m1_clump=m1_clump2, m1_freq=m1_freq)
    }
    return(environment())
}

###############################################################################
#'
#' Given regression parameters, compute a goodness-of-fit statistic of the
#' model with reg_par and data (tg_mat). This is useful for fitting
#' statistical models within optimization routines
#' 
#' The statistical model is 'largely' fit in 'numerical' space, not physical
#' space.
#' E.g. the smallest non-zero wavenumber resolvable by a DFT is 1/N,
#' and the nyquist wavenumber is 0.5
#' So be careful with units of kcx, kcy etc 
#'
#' @param reg_par = vector of 2 regression parameters  [kcx,kcy in numerical
#'        space]. A 3rd parameter may be accepted in some cases, see code
#' @param tg_rast = slip grid raster to fit
#' @param local_hurst = value for the hurst exponent
#' @param verbose = TRUE/FALSE -- Verbose error messages
#' @param default_seed= integer -- passed to set.seed for reproducible fitting
#'        with random fault generation (original .Random.seed is restored at the end)
#' @param NumRandSf = Number of slip distributions simulated to compute the
#'        goodness-of-fit of the model
#'
#' @return A goodness-of-fit measure -- minimising this will lead to the 'best'
#'         model fit
slip_goodness_of_fit<-function(
    reg_par,
    tg_rast,
    local_hurst = config_pars$hurst,
    verbose=FALSE,
    default_seed=1,   
    NumRandSf=200
    ){
    
    # Reproducible random seed 
    initial_seed = get_random_seed()
    # Reset the random seed to the original value when the function exits
    on.exit({.Random.seed <<- initial_seed})
    set.seed(default_seed) 

    #kcx=reg_par[1]
    #kcy=reg_par[2]
    #n=reg_par[3]

    # Treat invalid values -- we assume any reg_par[3] relates to 'alpha' in
    # the stable distribution
    if( (reg_par[1]<=0) | (reg_par[2]<=0)) return(9.0e+100)
    if( (length(reg_par)>2) && ((reg_par[3]<0)|reg_par[3] > 2)) return(9.0e+100)
   
    tg_mat = as.matrix(tg_rast)
    mean_tg_mat = mean(tg_mat)

    tmp = get_wavenumbers(tg_mat)
    kx = tmp[[1]]
    ky = tmp[[2]]


    # Store spectrum of all synthetic models here
    fake_f = array(NA,dim=c(nrow(kx), ncol(kx), NumRandSf))

    # Store fraction of fault covered with zeros here
    data_zero_frac = sum(tg_mat>0.0)/length(tg_mat)


    # Generate many random faults
    for(i in 1:NumRandSf){
        fake_data_cliprast = simulate_sffm(reg_par,tg_mat, 
            local_hurst = local_hurst)
        fake_f[,,i] = fake_data_cliprast
    }


    # Experiment with only keeping a zero-wavenumber component
    #trick_mat = tg_mat*0
    #trick_mat[,1] = 1
    #trick_mat[1,] = 1
    #trick_mat=1

    # Compute the 'goodness of fit' statistic
    spec_data = c(Mod(fft(tg_mat)))
    modstat = matrix(NA,ncol=length(spec_data),nrow=NumRandSf)
    for(i in 1:NumRandSf){
        modstat[i,] = c(Mod(fft(fake_f[,,i])))
    }
    mean_mod = colMeans(modstat)
    output_stat = mean((mean_mod-spec_data)**2)

    if(verbose){
        print(c(output_stat, reg_par, sum(output_stat), data_zero_frac))
    }

    return(output_stat)

}

#######################################################################
#'
#' Function to compute the statistical fit
#'
#' @param m1 matrix with slip values
#' @param local_hurst fixed non-default value for the hurst exponent
#' @param default_seed force the random seed value
#' @param NumRandSf number of synthetic slip simulations to use in goodness-of-fit computation
#' @return an object from 'optim' with the fit
#' 
fit_slip_parameters<-function(
    m1,
    local_hurst = config_pars$hurst,
    default_seed=1,   
    NumRandSf=400
    ){

    switch(config_pars$fourierFun,
        # 2 parameter model
        '2parSom' = {
            reg_par_start=c(2/dim(m1)[1], 2/dim(m1)[2])
        },
        # With variable hurst
        '3parSom' = {
            reg_par_start=c(2/dim(m1)[1], 2/dim(m1)[2], 1.0)
        },
        stop('config_pars$fourierFun not recognised')
    )
   

    # 3 parameter model
    #reg_par_start=c(2/dim(m1)[1], 2/dim(m1)[2],0.05)
    # 1 parameter model
    #reg_par_start=c(2/dim(m1)[1])
    
    # 3 parameter model
    #reg_par_start=c(2/dim(m1)[1], 2/dim(m1)[2], 0.4)
    m1_matrix = as.matrix(m1)
    output=try(
           optim(par=reg_par_start, 
                 fn=slip_goodness_of_fit,
                 tg_rast=m1_matrix,
                 local_hurst = local_hurst,
                 default_seed = default_seed,
                 NumRandSf = NumRandSf,
                 #hessian=TRUE,
                 control=list(maxit=2000))
            )

    if(class(output)=='try-error'){
        #print(output)
        output=list(convergence=-999, output=output, 
            reg_par_start=reg_par_start, m1=m1)

    }else if(output$convergence!=0){

        # Run again to make convergence more likely
        # (e.g. for Nelder/Mead, there might not have been enough iterations)
        reg_par_start=output$par
        output=optim(par=reg_par_start, 
                     fn=slip_goodness_of_fit,
                     tg_rast=m1_matrix,
                     local_hurst = local_hurst,
                     default_seed = default_seed,
                     NumRandSf = NumRandSf,
                     #hessian=TRUE,
                     control=list(maxit=2000))
    }

    return(output)
}

###############################################################################
#'
#' Extract the 'white noise' from a slip raster. It is not really white noise, since
#' it has constant amplitude spectrum (whereas a random image will have 
#' statistical variations in its amplitude spectrum). Also, it has constraints on
#' its variance and mean which mean the values are not iid, even if the original
#' image had iid values. Still, it is useful
#'
#' Note by construction if the mean slip is positive, this 'white-noise' will
#' have mean = 1/prod(dim(myrast)), and also (by Parsevals theorem)
#' sum(white_noise**2) = 1, so variance = 1/prod(dim(myrast)) assuming
#' mean(white_noise) is 0.  So it's not perfect white noise, but its
#' distribution does reflect what we should simulate to make the phase for SFFM
#' 
#' @param myrast raster or matrix
#' @return The inverse fourier transform of the exp(i*Arg(myrast)), where
#'         Arg(myrast) is the discrete fourier phase of myrast. 
#'
white_noise_extraction<-function(myrast){

    if(is.matrix(myrast)){
        ft = fft(myrast)
    }else{
        ft = fft(as.matrix(myrast))
    }

    # For R's inverse DFT, we must divide by the product of the image dim's
    white_noise = 1/prod(dim(ft))*fft(exp(1i*Arg(ft)), inverse=T)

    # Note -- if mean(myrast) is positive, the 'white_noise' above will always
    # have a mean = 1/prod(dim(ft)) (or negative of this if mean(myrast) is
    # negative). Remove that here.
    #white_noise = white_noise - mean(white_noise)

    stopifnot(all.equal(Im(white_noise), 0.*Im(white_noise)))
    return(Re(white_noise))

}


.test_white_noise_extraction<-function(){
    # Check that if we take an image and extract its 'white-noise',
    # then the argument of the extracted white noise is the same
    # as the argument of the original image     

    # Initial image [note: This will not have all spectra exactly = 1]
    wn = rnorm(100)
    dim(wn) = c(10, 10)

    wn2 = white_noise_extraction(wn)

    wn_Arg = Arg(fft(wn))

    stopifnot(all(abs(Arg(fft(wn2)) - wn_Arg)%%(2*pi) < 1.0e-12))

    print('PASS')

}
###############################################################################
#'
#' Wrapper code to compute the okada ground displacements
#'
#' The raster's do not give us enough information (e.g. we don't know dip/depth etc)
#' So we pass the name of the flt file used to generate them
#'
#' @param FFM_FLT filename of .flt or .FSP file with slip distribution
#' @param FFM_rast rasterLayer containing the rasterized slip surface from FFM_FLT
#' @param identical_rasters TRUE/FALSE. If true, FFM_FLT and FFM rast are
#'        assumed to refer to the same data. In that case, a test is done.
#'        Set to FALSE for simulated data
#' @param nxout number of x cells in the output raster
#' @param nyout number of y cells in the output raster
#' @param verbose -- print lots of things
#' @param strk see rake
#' @param dip see rake
#' @param rake if not contained in the FFM_FLT file, then these
#'        are used to define the bulk (planar) fault geometry
#' @param returnPtData return the point dataset defining the fault slip & geometry
#' @param return_xyDef write a file containing x,y, surface perturbation,
#'        with x,y coordinates being ALONG-STRIKE, ALONG SURFACE IN THE DOWN
#'        DIP DIRECTION
#' @param simID If return_xyDef is TRUE, append simID to the file name as an identifer
#' @return If returnPtData=FALSE, a Rasterlayer with ground deformation, 
#'         in ungeoreferenced projected  coordinates, with 0,0 = okada origin of 
#'         earthquake slip. If returnPtData=TRUE, a list with previous raster
#'         layer, and also the slip point data from the .flt or .FSP file
#'         
runOkada<-function(
    FFM_FLT,
    FFM_rast, 
    identical_rasters=TRUE,
    nxout=500,
    nyout=500,
    verbose=FALSE,
    strk=NA,
    dip=NA,
    rake=NA,
    returnPtData=FALSE, 
    return_xyDef=FALSE,
    simID='',
    runStamp=''){


    ###########################################################################

    # Read the FLT data (or FSP data)
    if(extension(FFM_FLT)=='.flt'){

        fltLines = readLines(FFM_FLT)
        gridInfoLine = grep('#Fault_', fltLines)
        if(length(gridInfoLine) != 1){
            stop('ERROR: FFM File format parsing error')
        }
        startDataLine = grep('#Lat', fltLines) 
        if(length(startDataLine)!=1){
            stop('ERROR: FFM File format parsing error')
        }

        # get Point data from flt file
        ptData1 = read.table(text=fltLines[startDataLine:length(fltLines)], 
            header=FALSE)

    }else if(extension(FFM_FLT)=='.FSP'){

        fspLines = readLines(FFM_FLT)
        fspLines = gsub('\t', '    ', fspLines) # Remove tabs
        startDataLine = grep("%    LAT      ",fspLines)+2
        if(length(startDataLine) != 1){
            stop('ERROR: FSP File format parsing error')
        }
        ptData1 = read.table(text=fspLines[startDataLine:length(fspLines)],
            header=FALSE)

        # Get strike and dip the 'easy-way'
        strike = rep(strk,length(ptData1[,1]))
        dip = rep(dip,length(ptData1[,1]))

        # Convert ptData1 to the same format as it would be with the flt file.
        # Lat Lon Depth Slip rake strike dip. Some files don't have rake,
        # careful
        if(length(ptData1[1,])>=7 & 
            length(grep('RAKE', fspLines[startDataLine]))>0){

            ptData1 = cbind(ptData1[,c(1:2,5:7)],strike,dip)

        }else{
            ptData1 = cbind(ptData1[,c(1:2,5:6)], rep(rake,length(ptData1[,1])),
                strike,dip)
        }
    
    }else{
        stop('Unknown ascii file type')

    }

    ###########################################################################

    # get Point data from raster
    ptData2 = rasterToPoints(FFM_rast)

    #@
    #@ Re-sort ptData2 so that it is aligned with ptData1
    #@
    # ptData1 and ptData2 are ordered differently. ptData1 is ordered with x
    # changing fast, then y changing more slowly.
    ptData2_resortTrick = ptData2[,1] + 100000*ptData2[,2]
    d2_to_d1 = sort(ptData2_resortTrick, index.return=T)$ix
    #d1_to_d2=sort(d2_to_d1,index.return=T)$ix
    ptData2=ptData2[d2_to_d1,]

    if(length(ptData1[,1]) != length(ptData2[,1])){
        stop('ERROR: Different number of points in FSP file vs raster')
    }

    # Check it has worked by comparing the slip values, if the rasters were
    # identical
    if(identical_rasters & 
       #(any(all.equal(ptData2[,3]*100, ptData1[,4])==FALSE))){
       (!isTRUE(all.equal(ptData2[,3], ptData1[,4])))){
         stop('ERROR in resorting finite fault data')
    }

    #@ Convert to x,y,z for okada function call.
    #@ 
    #@ Here we shift the coordinates to a projected system, with the
    #@ near-surface-centre-fault coordinate at the origin (i.e. top centre)
    #@
    
    nearSurf = which(ptData1[,3]==min(ptData1[,3])) # Indices of near surface points
    strk_C = ptData1[nearSurf[1],6]/180*pi # Take strike from strike of near-surface subfault
    # Find the corner point index with smallest along-strike coordinate
    tol = 1.0e-10

    if(strk_C < pi/2 -tol | strk_C > 3*pi/2+tol){
        # Corner has smallest latitude of all near surface points
        xcornerIndex = nearSurf[which.min(ptData1[nearSurf,1])]
    }else if(strk_C>pi/2+tol & strk_C < 3*pi/2-tol){
        # Corner has largest latitude of all near surface points
        xcornerIndex = nearSurf[which.max(ptData1[nearSurf,1])]
    }else if(all.equal(strk_C,pi/2,tolerance=tol)){
        # Corner has smallest longitude of near surface points
        xcornerIndex = nearSurf[which.min(ptData1[nearSurf,2])]
    }else if(all.equal(strk_C,3*pi/2,tolerance=tol)){
        # Corner has largest longitude of near surface points
        xcornerIndex = nearSurf[which.max(ptData1[nearSurf,2])]
    }else{
        stop('Cannot find xcornerIndex')
    }
        
    # Here is the coordinate system with x along strike, y down dip
    x = ptData2[,1] -(ptData2[xcornerIndex,1])
    # NOW SHIFT 'x' SO THE ORIGIN IS AT THE CENTRED AT THE TOP-CENTRE OF THE FAULT
    x = x-(res(FFM_rast)[1]*dim(FFM_rast)[2]/2)
    
    y_downdip = ptData2[,2]
    slip = ptData2[,3]
    # Get y [along surface] by Pythagoras
    y = sqrt( (y_downdip-min(y_downdip))**2 -(ptData1[,3]-min(ptData1[,3]))**2)

    # Rotate and flip coordinates to physical space 
    xunitRot = c(sin(strk_C), cos(strk_C))
    yunitRot = c(cos(strk_C), -sin(strk_C))

    newX = (x*xunitRot[1] + y*yunitRot[1])*1000.
    newY = (x*xunitRot[2] + y*yunitRot[2])*1000.

    #@ Various other parameters for okada
    strk = ptData1[,6] # degrees
    dip = ptData1[,7] # degrees
    lnth = rep(res(FFM_rast)[1],length(dip)) # km
    width = rep( res(FFM_rast)[2], length(dip)) # km
    rake = ptData1[,5]/180*pi 
    disl1 = slip*cos(rake) # m
    disl2 = slip*sin(rake) # m
    # Depth is computed for the cell centroid [offset from fsp file definition]
    depth = ptData1[,3] + width/2*cos(dip/180*pi) # km 

    # Check for no surface protrusions
    if(any(width/2*cos(dip/180*pi)>=depth)){
        stop('Earthquake sub-fault extends out of earth surface')
    }

    if(diff(range(width))>min(width)/10.){
        stop('Widths are not close to constant')
    }

    if(diff(range(lnth))>min(lnth)/10.){
        stop('lnths are not close to constant')
    }
  
    # Make a grid at which to get new raster values 
    # The surface deformation is computed on a rectangle this many times bigger
    # than the rupture, with smoothing at the extremities
    len_expansion_factor = 3.5     
    lscale = max(max(abs(newX)),max(abs(newY)))
    rlonVec = len_expansion_factor*seq(-lscale,lscale, len=nxout) # m
    rlatVec = len_expansion_factor*seq(-lscale,lscale, len=nyout) # m
    outgrid = expand.grid(rlonVec,rlatVec)

    if(verbose){
        print('')
        print('Some Coordinate extremes:')
        print(c(max(newX),min(newX),max(newY),min(newY)))
        print('Depth range:')
        print(c(min(depth),max(depth)))
        print('Strike, dip range :')
        print(c(min(strk), max(strk), min(dip), max(dip)))
        print('Outgrid bbox')
        print(c(min(outgrid[,1]), max(outgrid[,1]), min(outgrid[,2]), 
            max(outgrid[,2])))
        print('Disl1 range')
        print(c(min(disl1),max(disl1)))
        print('Disl2 range')
        print(c(min(disl2),max(disl2)))
        print('')
    }
    
    # generate tsunami with okada
    library(EqSim)
    ourTsunami = okada_tsunami(newX, newY, depth, strk, dip, lnth, width, 
        disl1, disl2, outgrid[,1],outgrid[,2], verbose=verbose)

    if(return_xyDef){
        #@
        #@ Write a file with the water surface perturbation, in along-strike /
        #@ along-surface-in-down-dip-direction coordinates with origin at top-centre

        # Rotate-flip x,y coordinates back to 'along-strike / perpendicular to
        # strike'
        # Note that the rotate-flip operation is its own inverse, so we can use
        # the old transformation
        newX2 = (outgrid[,1]*xunitRot[1] + outgrid[,2]*yunitRot[1])
        newY2 = (outgrid[,1]*xunitRot[2] + outgrid[,2]*yunitRot[2])

        #
        # Get xyz for deformation in ANUGA
        outxyz = cbind(newX2,newY2,ourTsunami[[3]])

        # To save file space, round deformation
        outxyz[,1] = round(outxyz[,1], 1)
        outxyz[,2] = round(outxyz[,2], 1)
        outxyz[,3] = round(outxyz[,3], 4)

        #@ Smooth the edges of the computed perturbation to zero
        #@ Otherwise, our tsunami deformation has a
        #@ 'seam' near the boundaries of the computed perturbation
        #@ Physically the 'seam' shouldn't matter, since the 'jump' is small
        #@ (mm) -- but it looks bad in ANUGA viewer
        
        outR = sqrt(outxyz[,1]**2+outxyz[,2]**2)
        outR = pmax(outR/(lscale*len_expansion_factor)*(1.25)-1,0.) # only > 0 near the edges of the computed perturbation
        outxyz[,3] = outxyz[,3]*pmax((0.25-outR)/0.25,0.)
      
        # Write outputs 
        scenario_base = gsub(extension(FFM_FLT), "", basename(FFM_FLT))
        dir.create(
            paste0(config_pars$deformation_folder,'/', runStamp, 
                '/', scenario_base), 
            recursive=TRUE, 
            showWarnings=FALSE)

        xyzFile = paste0(config_pars$deformation_folder, '/', runStamp, '/', 
            scenario_base, '/OceanInitial_', simID, '.xyz')
        write.table(outxyz,file=xyzFile, row.names=FALSE, col.names=FALSE, 
            sep=",")
    }

    # Return zdsp component as a raster
    outRast = rasterFromXYZ(cbind(outgrid,ourTsunami[[3]]),
        res=c(rlonVec[2]-rlonVec[1], rlatVec[2]-rlatVec[1]))

    if(returnPtData){
        return(list(outRast,ptData1)) 
    }else{
        return(outRast) 
    }
        
}

###############################################################################
#'
#' We make a plot of the ground surface deformation and compare
#' visually with a plot from a paper, which uses the same source
#'
#' The accuracy of the 'core' okada computation codes have been tested 
#' elsewhere -- this allows a check of the runOkada wrapping. 
#' Note the orientation, features, and magnuitude of
#' the deformation appear similar as in the paper (as far as we can tell).
#' 
check_runOkada<-function(){

    suppressPackageStartupMessages(library(raster))
    library(colorRampPC)
    blue2red=colorRampPC('cpt-city/kst/33_blue_red.cpt',n=40)
    
    comparison_file = paste0('../../DATA/SRCMOD_PROCESSED/',
        'SRCMOD_SUBDUCTION_DATA/slipRasts/',
        'S_2011003011_M9.1_Tohoku-Oki--Japan_13.tif')
    myrast=raster(comparison_file)
    myFSP=gsub('.tif','.FSP', comparison_file )
    mystrike=198
    mydip=10
    myrake=67.18

    deformation=runOkada(myFSP, myrast, strk=mystrike, dip=mydip, rake=myrake,
        verbose=TRUE, return_xyDef=TRUE, simID='Test')

    pdf(paste0(config_pars$quickfig_folder,'/GrilliComparison.pdf'),width=14,height=6)
    par(mfrow=c(1,2))
    par(oma=c(1,1,1,1))
    plot(myrast,xlab='Along Strike (km)', ylab='Down dip (km)',
         col=blue2red,asp=1)
    title('Shao (2011) source for Tohoku')
    plot(deformation[[1]], xlim=c(-337046.85, 78492.34),
         col=blue2red,
         ylim=  c(-305334.77, 305978.93), xlab='x',ylab='y',
         zlim=c(-maxValue(deformation[[1]]),maxValue(deformation[[1]])), 
          asp=1)
    title('Compare with Fig 3 in Grilli et al., 2012')
    par(mfrow=c(1,1)) 
    grillifig=brick('../../DATA/TEST_DATA/Grilli_Fig3.png')
    plotRGB(grillifig)
    dev.off()
}

##############################################################################
#'
#' Plot the ocean deformation xyDef file results
#'
#' @param xyDef is a 3 column matrix with x,y, ocean surface deformation
#' @param ... other parameters passed to 'plot'
#'
plot_xyDef<-function(xyDef,...){
    library(colorRampPC)

    pos_col=colorRampPC('cpt-city/nd/basic/Cyan_Magenta_CW.cpt',n=100)
    #pos_col=colorRampPC('cpt-city/kst/22_hue_sat_value2.cpt' ,n=100)
    #pos_col=colorRampPC('cpt-city/kst/33_blue_red.cpt',n=100)
    #pos_col=colorRampPC('cpt-city/nd/basic/Blue_Cyan_CCW.cpt',n=100) 
    neg_col=colorRampPC('cpt-city/nd/atmospheric/Sunset_Real.cpt',n=100)

    mycol=xyDef[,3]*NA
    mycol[xyDef[,3]<0] = neg_col[as.numeric(cut(xyDef[xyDef[,3]<0,3], 100))]
    mycol[xyDef[,3]>=0] = pos_col[as.numeric(cut(xyDef[xyDef[,3]>=0,3], 100))]

    par(oma=c(0,0,0,2))
    layout(matrix(c(1,2),ncol=2),widths=c(0.8,0.2))
    par(mar=c(4,4,2,0))
    plot(xyDef[,1:2],pch='.',col=mycol,asp=1,...)

    myclz=sort(xyDef[,3],index.return=T)
  
    par(mar=c(1,0,0,0))  
    newColVec=resample_colorVec(mycol[myclz$ix],
        breaks=c(myclz$x,max(myclz$x)+0.001), n=400)
    newbrks=seq(min(myclz$x),max(myclz$x)+0.001,len=401)
    plot_colorVec(newColVec, breaks=newbrks,vertical=TRUE,
        add_axis=TRUE,plotWidthScale=0.3,plotHeightScale=0.5)

}

###############################################################################
#'
#' Plot the initial deformation xyDef files in each folder
#'
#' @param defFolder a folder with lots of .xyz files containing the deformation
#' @param ... other parameters passed to plot_xyDef
#'
plot_def_folder<-function(defFolder, ...){
    inFiles=dir(defFolder,pattern='.xyz',full.names=T)
    if(length(inFiles)==0) return
    pdf(paste0(dirname(inFiles[1]),'/plots.pdf'),width=8,height=5)
    for(ff in inFiles){
        inDat=read.table(ff,sep=",")
        plot_xyDef(inDat,main=basename(ff), cex=2, ...)
    }
    dev.off()
    return
}

#' Make a quick plot of ocean deformation, similar to plot_xyDef
#'
#' @param okadaDef result of runOkada
#' 
plot_okadaDef<-function(okadaDef){

    library(colorRampPC)

    pos_col=colorRampPC('cpt-city/nd/basic/Cyan_Magenta_CW.cpt',n=100)
    neg_col=colorRampPC('cpt-city/nd/atmospheric/Sunset_Real.cpt',n=100)

    MM=maxValue(okadaDef)
    mm=minValue(okadaDef)

    colVal=c(neg_col,pos_col)
    colBreaks=c(seq(mm,-1e-06,len=100),seq(0,MM,len=100),MM+1.0e-06)

    oldMar=par('mar')
    newMar=c(oldMar[1:3],0)
    par('mar'=newMar)
    layout(matrix(c(1,2),ncol=2),widths=c(0.75,0.25))
    image(okadaDef, col=colVal, breaks=colBreaks, asp=1)

    newCol=resample_colorVec(colVal,colBreaks,n=200)
    newColBreaks=seq(min(colBreaks),max(colBreaks),len=201)
    par('mar'=c(oldMar[1],0,oldMar[3:4]))
    plot_colorVec(newCol,breaks=newColBreaks,add_axis=TRUE,vertical=TRUE)
   
    par('mar'=oldMar) 
    
}

