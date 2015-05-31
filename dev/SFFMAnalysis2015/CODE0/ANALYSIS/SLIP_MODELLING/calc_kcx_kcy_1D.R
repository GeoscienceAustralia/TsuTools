
calc_1D_kcx_kcy<-function(tg,fitmethod='zero_wavenumbers', PLOT=FALSE){
    ##
    ## Fit the '1-dimensional' estimates of kcx and kcy
    ## If the model has been clipped / abs'ed / tapered etc, this approach can be bias,
    ## used here for comparison purposes only
    ##
    ## INPUT: tg = raster of finite fault
    ##        fitmethod = Method to use, see the code  
    ##
    ## WARNING: The outputs are 'NUMERICAL' wavenumbers. To get to physical wavenumbers,
    ##       you have to divide kcx (kcy) by the x (y) raster cell size
    ##
    ##

    tg_mat = as.matrix(tg)
    tg_fft_spec = Mod(fft(tg_mat))

    # Make matrices with x,y wavenumbers
    B1 = nrow(tg_fft_spec)
    A1 = ncol(tg_fft_spec)
    X2 = pmin(0:(B1-1), B1-(0:(B1-1)))/B1
    X1 = pmin(0:(A1-1), A1-(0:(A1-1)))/A1
    kx = matrix(X1,ncol=A1,nrow=B1, byrow=T)
    ky = matrix(X2, ncol=A1,nrow=B1)

    switch(fitmethod, 
    '1D_slices'={
        ## Compute 1D fourier spectra in x (and then y) directions, average,
        ## and fit to that
        ## This doesn't seem to work
        x_fft_spec_1d = NA*tg_fft_spec
        y_fft_spec_1d = NA*tg_fft_spec
        for(i in 1:B1){
            x_fft_spec_1d[i,] = Mod(fft(tg_mat[i,]))
        }
        for(i in 1:A1){
            y_fft_spec_1d[,i] = Mod(fft(tg_mat[,i]))
        }

        x_mean_spec = colSums(x_fft_spec_1d)/B1
        y_mean_spec = rowSums(y_fft_spec_1d)/A1
    }, 
    'zero_wavenumbers'={
        # This should work if the fourier spectrum agrees with the model
        x_mean_spec = tg_fft_spec[1,]
        y_mean_spec = tg_fft_spec[,1]

        # Add a perturbation (to allow the nonlinear fitting alg to converge)
        x_mean_spec = x_mean_spec+1e-2*min(x_mean_spec)*(runif(length(x_mean_spec))-0.5)
        y_mean_spec = y_mean_spec+1e-2*min(y_mean_spec)*(runif(length(y_mean_spec))-0.5)
    }

    )

    kx1D = kx[1,]
    ky1D = ky[,1]


    ## Fit x
    x_mean_spec_No_0 = (x_mean_spec/max(x_mean_spec)) 
    kx1D_No_0 = kx1D 
    Y = x_mean_spec_No_0
    X = kx1D_No_0
    modx = nls(log10(Y)~log10((1+(X/alpha)**4)**(-0.5)),
               start=list(alpha=1/A1))

    alpha = coef(modx)[1]
    beta = 1.0 
    kcx = c(alpha)
   
    ## Plot 
    if(PLOT){
        par(mfrow=c(2,1))
        plot(kx1D, x_mean_spec_No_0,ylim=c(min(x_mean_spec_No_0)*0.3, max(x_mean_spec_No_0)),log='y')
        kx1Dcon = seq(0,max(kx1D),len=100)
        points(kx1Dcon, beta*(1+(kx1Dcon/alpha)**4)**(-0.5),
               t='l',col='red')
    }

    ## Fit y
    y_mean_spec_No_0 = (y_mean_spec/max(y_mean_spec))
    ky1D_No_0 = ky1D
    Y = y_mean_spec_No_0
    X = ky1D_No_0
    mody = nls(log10(Y)~log10((1+(X/alpha)**4)**(-0.5)),
               start=list(alpha=1/B1))
    alpha = coef(mody)[1]
    beta = 1.0 

    ## Plot
    if(PLOT){
        plot(ky1D, y_mean_spec_No_0,ylim=c(0.3*min(y_mean_spec_No_0), max(y_mean_spec_No_0)))
        ky1Dcon=seq(0,max(ky1D),len=100)
        points(ky1Dcon, beta*(1+(ky1Dcon/alpha)**4)**(-0.5),
               t='l',col='red')
    }
    kcy = c(alpha)

    out = c(kcx,kcy)
    names(out) = c('kcx', 'kcy')
    return(out)

}

