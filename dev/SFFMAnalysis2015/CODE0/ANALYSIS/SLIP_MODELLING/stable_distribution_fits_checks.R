# Code to test / fit stable distributions
library(raster)
library(parallel)
# Ensure appropriate parallel random number generator is used (or the numbers
# might not be independent)
RNGkind("L'Ecuyer-CMRG")


# Make fake configuartion parameters to keep stable.R happy
if(!exists('config_pars')){
    config_pars = list(
        MC_CORES=12, 
        myrastlist = Sys.glob( '../../DATA/SRCMOD_PROCESSED/SRCMOD_SUBDUCTION_DATA/slipRasts/*.tif' ),
        pdf_folder = 'OUTPUTS_2/pdfs',
        rimages_folder = 'OUTPUTS_2/Rimages'
        )
       
    config_pars$myrasts = lapply( config_pars$myrastlist, f<-function(x) raster(x))

    # Extract Mw from the filename (hacky)
    config_pars$mw = as.numeric(
        gsub('M', '', 
            gsub('_','.', 
                unlist(
                    lapply(strsplit(basename(config_pars$myrastlist), '_'), 
                        f<-function(x) x[3])
                )
             )
        )
    )
}

source('stable.R')

if(TRUE){

    # This code repeats the fit test 50 times
    fitting_tests = mclapply(as.list(1:50), test_stable_dist_fit, 
        mc.cores=config_pars$MC_CORES, mc.preschedule=FALSE)

    #############################################################
    # Check the cauchy fit parameter estimates 
    ############################################################
    cauchy_fit_par = matrix(unlist(lapply(fitting_tests, 
        f<-function(x){
            out = x$cauchy_fit$estimate
            if(length(out)!=4) out = rep(NA, 4)
            return(out)}
             )), ncol=4, byrow=T)
    # Plot it
    par(mfrow=c(2,2))
    true_par = c(1, 0, 1, 0)
    for(i in 1:4){
        hist(cauchy_fit_par[,i])
        abline(v = true_par[i], col='red')
    }

    # Stats on confidence intervals and failure to converge
    # Note these confidence intervals are asymptotic only, might not be very good
    # See instead the profile likelihood CI
    cauchy_fit_check = mean(unlist(lapply(fitting_tests, 
        f<-function(x) all(abs(x$cauchy_fit$estimate - true_par)/(2*(x$cauchy_fit$sd)) < 1.0) ))
        , na.rm=T)
    cauchy_fit_failure = mean(is.na(unlist(lapply(fitting_tests, 
        f<-function(x) all(abs(x$cauchy_fit$estimate - true_par)/(2*(x$cauchy_fit$sd)) < 1.0) ))))

    #############################################################
    # Check the stable_1 parameter estimates
    #############################################################
    stable_1_fit_par = matrix(unlist(lapply(fitting_tests, 
        f<-function(x){
            out = x$stable_1_fit$estimate
            if(length(out)!=4) out = rep(NA, 4)
            return(out)}
             )), ncol=4, byrow=T)

    par(mfrow=c(2,2))
    true_par = c(1.7, 0.4, 10, 0)
    for(i in 1:4){
        hist(stable_1_fit_par[,i])
        abline(v = true_par[i], col='red')
    }

    # Stats on confidence intervals and failure to converge
    stable_1_fit_check = mean(unlist(lapply(fitting_tests, 
        f<-function(x) all(abs(x$stable_1_fit$estimate - true_par)/(2*(x$stable_1_fit$sd)) < 1.0) ))
        , na.rm=T)
    stable_1_fit_failure = mean(is.na(unlist(lapply(fitting_tests, 
        f<-function(x) all(abs(x$stable_1_fit$estimate - true_par)/(2*(x$stable_1_fit$sd)) < 1.0) ))))

    ############################################################
    # Check the stable_2 parameter estimates
    ############################################################
    stable_2_fit_par = matrix(unlist(lapply(fitting_tests, 
        f<-function(x){
            out = x$stable_2_fit$estimate
            if(length(out)!=4) out = rep(NA, 4)
            return(out)}
             )), ncol=4, byrow=T)
    par(mfrow=c(2,2))
    true_par = c(0.8, -0.1, 1, -10)
    for(i in 1:4){
        hist(stable_2_fit_par[,i])
        abline(v = true_par[i], col='red')
    }
    # Stats on confidence intervals and failure to converge
    stable_2_fit_check = mean(unlist(lapply(fitting_tests, 
        f<-function(x) all(abs(x$stable_2_fit$estimate - true_par)/(2*(x$stable_2_fit$sd)) < 1.0) ))
        , na.rm=T)
    stable_2_fit_failure = mean(is.na(unlist(lapply(fitting_tests, 
        f<-function(x) all(abs(x$stable_2_fit$estimate - true_par)/(2*(x$stable_2_fit$sd)) < 1.0) ))))


}

# Here we fit the stable distribution to the FFI
# NOTE: Our SFFM can potentially have phases with a distribution that
#       differs from the distribution used to generate the original 
#       phase generating field (because of post-processing). So rather 
#       than using direct fits, we choose alpha
fits_to_FFI = fit_all_rasters_with_stable_dist_white_noise(config_pars$myrasts,
    symmetric=TRUE)
plot_fitted_stable_dist_white_noise(fits_to_FFI, config_pars$myrasts)
stable_pars = get_stable_dist_all_parameter_estimates(fits_to_FFI)

# Print the various parameters + do a basic plot
summary(stable_pars[[1]])
plot_stable_dist_parameter_estimates(fits_to_FFI, config_pars$mw)

# Get better confidence intervals than the default fitdist inverse-hessian ones
profile_fits_to_FFI = mclapply(fits_to_FFI, profile_stablesym, 
    make_plot=FALSE, mc.preschedule=FALSE, mc.cores=config_pars$MC_CORES)

fitted_alpha_vals = unlist(lapply(profile_fits_to_FFI, f<-function(x) x$best_fit_pars[1]))
median_fitted_alpha = median(fitted_alpha_vals)
# Coverage in profile
median_coverage = unlist(lapply(profile_fits_to_FFI, f<-function(x) (prod(x$alpha_range - median_fitted_alpha) < 0.) ) )
# %age coverage
mean(median_coverage)

save.image(paste0(config_pars$rimages_folder, '/stable_distribution_fits.Rdata'))
