###############################################################################
# Routines for fitting / testing fits of stable distributions in our context
#
# AUTHORS
# Gareth Davies, Geoscience Australia, 2013-2015,
# gareth.davies.ga.code@gmail.com
#
###############################################################################


library(stabledist)
library(fitdistrplus)

if(!exists('white_noise_extraction', mode='function')){
    # Get the 'white_noise_extraction' function (which is not really white!)
    source('fit_simulate_earthquake.R') 
}

###############################################################################
#'
#' Wrapper to stabledist::dstable -- to prevent failure with invalid parameters
#' Only required for dstable since that is used in ML optimization.
#' See ?stabledist::dstable for parameters
dstable<-function(x, alpha, beta, gamma = 1, delta = 0, pm = 0,
    log = FALSE, tol = 64*.Machine$double.eps, zeta.tol = NULL,
    subdivisions = 1000){

    output = try(
        stabledist::dstable(x, alpha, beta, gamma, delta, pm, log, tol, 
            zeta.tol, subdivisions),
        silent=TRUE
        )

    # Here we hack around the errors for invalid parameters in
    # stabledist::dstable
    if(class(output)=='try-error'){
        if(log == FALSE){
            # Very small number, but the log transform is not -Inf
            return(rep(1.0e-300, length(x)))
        }else{
            return(rep(log(1.0e-300), length(x)))
        }
    }

    return(output)

}


#' Negative log likelihood for stable distribution
#' See dstable for more info
#'
#' @param pars = c(alpha, beta, gamma, delta0)
negloglik_dstable<-function(pars, data){
    output = -sum(dstable(data, alpha=pars[1], beta=pars[2], gamma=pars[3], 
        delta=pars[4], log=TRUE))

    return(output)
}

###############################################################################
#
# d,p,q,r functions for the SYMMETRIC stable distribution (beta=0,delta=0)


#' Density function for the symmetric stable distribution (beta=0, delta=0)
#' See ?stabledist::dstable for parameters
dstablesymm<-function(x, alpha, gamma=1, pm = 0,
    log = FALSE, tol = 64*.Machine$double.eps, zeta.tol = NULL,
    subdivisions = 1000){

    output = dstable(x, alpha, beta=0, gamma, delta=0, pm,
        log, tol, zeta.tol, subdivisions)
    
    return(output)
}

#' CDF for the symmetric stable distribution (beta=0, delta=0)
#' See ?stabledist::pstable for parameters
pstablesymm<-function(q, alpha, gamma = 1, pm = 0,
             lower.tail = TRUE, log.p = FALSE,
             tol = 64*.Machine$double.eps, subdivisions = 1000){
    
    #pstable(q, alpha, beta, gamma = 1, delta = 0, pm = 0,
    #         lower.tail = TRUE, log.p = FALSE,
    #         tol = 64*.Machine$double.eps, subdivisions = 1000)

    output = pstable(q, alpha, beta = 0, gamma, delta = 0, pm,
        lower.tail, log.p, tol, subdivisions)

    return(output)

}

#' Quantile function for the symmetric stable distribution (beta=0, delta=0)
#' See ?stabledist::qstable for parameters
qstablesymm<-function(p, alpha, gamma = 1, pm = 0,
             lower.tail = TRUE, log.p = FALSE,
             tol = .Machine$double.eps^0.25, maxiter=1000,
             integ.tol = 1e-7, subdivisions = 200 ){

     #qstable(p, alpha, beta, gamma = 1, delta = 0, pm = 0,
     #        lower.tail = TRUE, log.p = FALSE,
     #        tol = .Machine$double.eps^0.25, maxiter = 1000,
     #        integ.tol = 1e-7, subdivisions = 200)

    output = qstable(p, alpha, beta = 0, gamma, delta = 0, pm,
        lower.tail, log.p, tol, maxiter, integ.tol, subdivisions)

    return(output)

}

#' Random sample function for symmetric stable distribution (beta=0, delta=0)
#' See ?stabledist::rstable for parameters
rstablesymm<-function(n, alpha, gamma=1, pm=0){
    
     # rstable(n, alpha, beta, gamma = 1, delta = 0, pm = 0)
    output = rstable(n, alpha, beta=0, gamma, delta=0, pm)
    return(output)
}

#' Negative log likelihood for symmetric stable distribution
#' See dstablesymm for more info
#'
#' @param pars = c(alpha, gamma)
negloglik_dstablesymm<-function(pars, data){

    output = -sum(dstable(data, alpha=pars[1], beta=0., gamma=pars[2], 
        delta=0., log=TRUE))

    return(output)
}

###############################################################################
#'
#' Profile around a best-fit to a symmetric stable distribution to find a
#' confidence interval for alpha
#'
#' Canned profiling methods didn't seem to work for these fits
#'
#' @param fit a maximum likelihood fit from fitdistr
#' @param ci_level level for profile confidence region (qchisq(ci_level,1)/2)
#' @param make_plot make a plot
#' @param nr, nc number of rows/columns for grid
profile_stablesym<-function(fit, ci_level=0.95, 
    make_plot=TRUE, nr=21, nc=20, verbose=TRUE){

    best_fit_pars = fit$estimate
    data = fit$data

    # Apparent 'best' negative log likelihood
    min_nll = negloglik_dstablesymm(best_fit_pars, data=data)

    # The bound of the negative log likelihood in the confidence interval
    desired_nll = min_nll + qchisq(ci_level, 1)/2

    # Compute nll on a grid
    store_nll = matrix(NA, nrow = nr, ncol = nc)
    alphas=seq(0.5*best_fit_pars[1], min(2.0*best_fit_pars[1], 1.9999), len=nr)
    gammas=seq(0.5*best_fit_pars[2], 1.5*best_fit_pars[2], len=nc)
    for(i in 1:nr){
        for(j in 1:nc){
            store_nll[i,j] = negloglik_dstablesymm(c(alphas[i], gammas[j]), 
                data=data)
        }
        if(verbose) print(i)
    }

    # Compute a contour at the profile confidence level
    profileContour = contourLines( alphas, gammas, store_nll, 
        levels = desired_nll)

    alpha_range = range(profileContour[[1]]$x)
    gamma_range = range(profileContour[[1]]$y)

    if(make_plot){
        image(alphas, gammas, store_nll, col=rainbow(20))
        contour(alphas, gammas, store_nll, levels = desired_nll, add=T)
    }
   
    return(list(alpha=alphas, gamma=gammas, nll = store_nll, min_nll = min_nll, 
        desired_nll = desired_nll, profileContour = profileContour, 
        alpha_range = alpha_range, gamma_range = gamma_range, 
        ci_level = ci_level, best_fit_pars = best_fit_pars, data = data)) 

}

###############################################################################
#'
#' Wrapper to fit the stable distribution, with default starting parameters
#' which are appropriate for the current application, but might not work
#' elsewhere
fit_stable_dist<-function(
    x, 
    symmetric = FALSE,
    start=list(alpha=1.5  , beta=0     ,gamma=5, delta=0),
    apply_limits=TRUE,
    optim.method='Nelder-Mead'){
    # Wrapper to fit a stable distribution
    # See warnings about accuracy of stabledist::dstable

    if(symmetric){
        # 2 parameters, alpha and gamma, since beta = 0 and delta = 0 for a
        # symmetric distribution
        dist = 'stablesymm'
        start = list(alpha = start$alpha, gamma=start$gamma)
        # Theoretical ranges:
        # alpha [0,2]  (but numerics might be bad near alpha=0)
        # gamma [0, Inf]
        lower = c(   0.01     ,   0)
        upper = c(2.0-1.0e-05 , Inf)
    }else{
        dist = 'stable'
        # Theoretical ranges:
        # alpha [0,2]   (but numerics might be bad near alpha=0) 
        # beta [-1, 1]
        # gamma [0, Inf]
        # delta [-Inf, Inf]
        lower = c(   0.01     ,-1.0+1.0e-05,   0   , -Inf)
        upper = c(2.0-1.0e-05 , 1.0-1.0e-05, Inf   ,  Inf)
    }

    if(apply_limits){
        # Use L-BFGS to fit within a box (cannot have NA values)
        fit = fitdist(x, dist, start= start, lower=lower, upper=upper)

    }else{
        # No constraints (but we have tweaked dstable/dstablesymm to give
        # very small values outside the allowed range
        fit = fitdist(x, dist, start=start, optim.method=optim.method)

    }

    return(fit)
}

#' Wrapper to fit stable distribution to white noise from the phase of myrast
fit_stable_dist_to_white_noise<-function(myrast, symmetric=FALSE){

    white_noise = white_noise_extraction(myrast)

    fit_output = fit_stable_dist(
        c(white_noise), 
        symmetric = symmetric,
        apply_limits = FALSE, 
        optim.method='BFGS')

    return(fit_output)
}

#' 
#' Some checks that our stable distribution fit is working
#' This can be run many times (e.g in parallel) to check that the results are ok
#' @param i Dummy argument (not required to run, but mclapply needs an argument)
test_stable_dist_fit<-function(i=NULL){

    # Cauchy distribution -- alpha = 1, beta = 0
    cauchy_data = rcauchy(100)
    cauchy_fit = fit_stable_dist(cauchy_data, apply_limits = FALSE, 
        optim.method='BFGS') 

    # Stable distribution test 1
    stable_data_1 = rstable(100, alpha=1.7, beta = 0.4, gamma = 10, delta=0)
    stable_1_fit = fit_stable_dist(stable_data_1, apply_limits = FALSE, 
        optim.method='BFGS') 


    # Stable distribution test 1
    stable_data_2 = rstable(100, alpha=0.8, beta = -0.1, 
        gamma = 1, delta = -10)
    stable_2_fit = fit_stable_dist(stable_data_2, apply_limits = FALSE, 
        optim.method='BFGS') 

    # FFI of to
    tokachi_oki_FFI = '../../DATA/SRCMOD_PROCESSED/SRCMOD_SUBDUCTION_DATA/slipRasts/S_2003009025_M8.21_Tokachi-oki-Japan_51.tif'
    FFI = raster(tokachi_oki_FFI)
    stable_FFI_fit = fit_stable_dist_to_white_noise(FFI)
   
    return(environment()) 
}

#' Fit white noise from all rasters with a stable distribution, in parallel
fit_all_rasters_with_stable_dist_white_noise<-function(myrasts, symmetric=FALSE){
    library(parallel)
    white_noise_fits = mclapply(myrasts, fit_stable_dist_to_white_noise, 
        symmetric=symmetric, mc.cores = config_pars$MC_CORES, 
        mc.preschedule=FALSE)
    return(white_noise_fits)
}

#' Plot outputs of fit_all_rasters_with_stable_dist_white_noise
plot_fitted_stable_dist_white_noise<-function(white_noise_fits, myrasts){

    pdf(paste0(config_pars$pdf_folder,'/white_noise_stable_distribution.pdf'),
        width=10, height=15)

    for(i in 1:length(white_noise_fits)){

        fit = white_noise_fits[[i]]
        if(class(fit) == 'try-error') next

        # plot
        par(mfrow=c(4,2))
        plot(myrasts[[i]])
        title(myrasts[[i]]@file@name)
        image(white_noise_extraction(myrasts[[i]]), main='"white"-noise')
        # Normal qq plot
        qqnorm(fit$data)
        qqline(fit$data)

        # Diagnostics for stable distribution fit
        ppcomp(fit, main='P-P Plot (stable dist)')

        # the quantile function can struggle sometimes
        #try(qqcomp(fit, main='Q-Q Plot (stable dist)'))
        fitted_quantiles = sort(fit$data)
        prob_pts = ppoints(length(fitted_quantiles))
        # Can fail if prob_pts == 0.5 ! FIXME: Report as bug
        adjust = which(abs(prob_pts-0.5) < 1.0e-06)
        if(length(adjust)>0){
            prob_pts[adjust] = prob_pts[adjust] + 1.0e-03
        }
        if(length(fit$estimate)==2){
            theoretical_quantiles =try( qstablesymm(prob_pts, alpha=fit$estimate[1],
                gamma=fit$estimate[2]))
            if( class(theoretical_quantiles) == 'try-error'){
                browser()
            }
        }
        plot(theoretical_quantiles, fitted_quantiles)
        abline(0, 1, col='red')

        cdfcomp(fit, main='CDF Plot (stable dist)')
        denscomp(fit, main='Density Plot (stable dist)')
        
        plot(c(0,1),c(0,1),col=0,main=fit$estimate[1], cex.main=3)
    }    

    dev.off()
}

#' Extract the parameter estimates + their se from the white_noise_fits
get_stable_dist_all_parameter_estimates<-function(white_noise_fits){

    stable_estimates = lapply(white_noise_fits, 
        f<-function(x){
            if(class(x)!='try-error'){
                return(x$estimate) 
            }else{
                return(NA)
            }
        }
        )

    stable_sd = lapply(white_noise_fits, 
        f<-function(x){
            if(class(x)!='try-error'){
                return(x$sd) 
            }else{
                return(NA)
            }
        }
        )

    # Make sure NA values have the right length
    l=0
    for(i in 1:length(stable_estimates)){
        l = max(l, length(stable_estimates[[i]]))
    }


    for(i in 1:length(stable_estimates)){
        if( (length(stable_estimates[[i]]) == 1) && is.na(stable_estimates[[i]]) ){
            stable_estimates[[i]] = rep(NA, l)
        }
        if( (length(stable_sd[[i]])==1) && is.na(stable_sd[[i]]) ){
            stable_sd[[i]] = rep(NA,l)
        }
    }


    # Reformat
    stable_estimates = matrix(unlist(stable_estimates), ncol = l, byrow=T)
    colnames(stable_estimates) = names(white_noise_fits[[1]]$estimate)
    stable_sd = matrix(unlist(stable_sd), ncol = l, byrow=T)
    colnames(stable_sd) = names(white_noise_fits[[1]]$sd)

    return(list(stable_estimates, stable_sd))
}

#' Quick plot of parameter estimates + standard deviation
#' Optionally plot against provided Mw
plot_stable_dist_parameter_estimates<-function(white_noise_fits, mw = NULL){

    parameter_fits = get_stable_dist_all_parameter_estimates(white_noise_fits)

    par(mfrow=c(2,2))
    npar = length(parameter_fits[[1]][,1])
    parnames = colnames(parameter_fits[[1]])
   
    for(i in 1:length(parnames)){ 
        # Plot CI bars with 2 sd -- note this is only asymptotically correct,
        # doesn't consider boundaries
        if(is.null(mw)){
            xvals = 1:npar
        }else{
            xvals = mw
        }

        lower_ci = parameter_fits[[1]][,i]-2*parameter_fits[[2]][,i]
        upper_ci = parameter_fits[[1]][,i]+2*parameter_fits[[2]][,i] 
        plot(0, 0, xlim=range(xvals), 
            ylim=c(min(lower_ci, na.rm=T), max(upper_ci, na.rm=T)), col=0)
        arrows(xvals, lower_ci,  
            xvals, upper_ci,
            angle = 90, 
            code = 3)
        points(xvals, parameter_fits[[1]][,i])
        title(paste0(parnames[i], ' \n[asymptotic CI only, use profile likelihood instead]'))
        abline(h = median(parameter_fits[[1]][,i], na.rm=T), col='red', 
            lty='dashed')
    }
}
