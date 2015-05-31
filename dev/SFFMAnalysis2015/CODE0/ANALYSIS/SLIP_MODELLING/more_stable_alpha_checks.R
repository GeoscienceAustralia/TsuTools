# Here we check values of alpha

library(parallel)
library(methods)


#' Function to compute a quantile confidence interval
#'
#' @param x vector of numeric data
#' @param Quant the quantile we want a CI for
#' @param level the CI level
#'
quantile_CI<-function(x, Quant=0.5, level=0.95){

    alpha = (1-level)
    n = length(x) # Size of data

    allQuants = pbinom(0:n, n, Quant) # All quantiles
    lower_index = max(which(allQuants<=alpha/2))
    upper_index = min(which(allQuants>=(1-alpha/2)))

    if(lower_index==-Inf | upper_index==n+1){
        message = paste("Cannot compute Quantile ",
            Quant, " at " , level*100, " % confidence " ,
            " with " , n, " data points " )
        stop(message)
    }

    output=list()

    # Interval
    output$ci=sort(x)[c(lower_index,upper_index)]

    # Exact probability (>level)
    output$p = allQuants[upper_index]-allQuants[lower_index]

    return(output)
}


# Identify stable distribution model runs
stable_dist_runs = Sys.glob('OUTPUTS_2/Rimages/*_S_S*')

# Make place for text outputs
dir.create('alpha_output_results', showWarnings=FALSE)

# for(stable_dist_run in stable_dist_runs){

# Make a function to do our computation
test_alpha_fit_and_fit_to_ffi<-function(stable_dist_run){

    # Capture the output
    sink(paste0('alpha_output_results/', basename(stable_dist_run), '.txt'))

    print('Loading session')
    # We need to load to .GlobalEnv for everything to work
    load(stable_dist_run, envir=.GlobalEnv)

    print('Loading key files')
    # Not sure whether we need local=TRUE here or not, but it is working
    source('stable.R', local=TRUE)
    source('fit_simulate_earthquake.R', local=TRUE)



    ###########################################################################
    # For a range of alpha values, fit kcx/kcy for each SFFM. 

    # alpha in 1, 1.1, 1.2, ..., 1.9, (2.0 - epsilon)
    #alphas = c(seq(1.0, 1.9, len=10), 1.9999)

    #fit_pars_for_varying_alpha<-function(i, verbose=FALSE){

    #    if(verbose) print(i)

    #    real_rast_mat = as.matrix(rod$myrasts[[i]])

    #    # For every alpha, compute the best fit parameters for the data
    #    synthetic_mats = list()
    #    best_pars = list()
    #    for(k in 1:length(alphas)){
    #        # For alpha[k], find the best kcxN, kcyN
    #        f_to_optim<-function(best_k){
    #            slip_goodness_of_fit(
    #                c(best_k, alphas[k]), tg_rast = real_rast_mat, 
    #                local_hurst=1.0, default_seed = 1, NumRandSf = 400)
    #        }

    #        best_pars[[k]] = optim(c(0.1, 0.1), f_to_optim)

    #        if(verbose) print(best_pars[[k]])

    #    }

    #    return(environment())
    #}
    #
    #parallel_fits = mclapply(1:length(rod$myrasts), fit_pars_for_varying_alpha, 
    #    mc.preschedule=FALSE, mc.cores=12)

    #results_storage[[stable_dist_run]] = parallel_fits


    #' Compute the kurtosis of x
    #' (There are a range of definitions of this)
    #'
    kurtosis<-function(x){
        n = length(x)
        n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
    }

    #' Compute a goodness of fit measure for the stable distribution alpha parameter,
    #' assuming slip_mat was generated from our SFFM algorithm with that alpha parameter
    #' 
    #' @param alpha the alpha to get goodness of fit for
    #' @param slip_mat the slip matrix to fit
    #' @param wavenumber_par the corner wavenumber parameters assumed used to
    #'        generate slip_mat. These are often not very sensitive to alpha so can
    #'        be estimated by separate fitting
    #' @param nreps Replications in the stochastic algorithm used to generate teh goodness of fit measure
    #' @param default seed Deterministic random seed used in the stochastic goodness of fit algorithm
    #' @return The goodness of fit measure
    #'
    alpha_goodness_of_fit<-function(alpha, slip_mat, wavenumber_par, 
        nreps = 50, default_seed = 1234){
        
        wn_data = white_noise_extraction(slip_mat)

        # Use deterministic random numbers
        # Reset the random seed to the original value when the function exits
        initial_seed = get_random_seed()
        on.exit({.Random.seed <<- initial_seed})
        set.seed(default_seed) 

        stat = rep(NA, nreps)
        for(i in 1:nreps){
            # Simulate a model, extract its white noise, compute the kurtosis
            wn_random = simulate_sffm(c(wavenumber_par, alpha), slip_mat)
            wn_random = white_noise_extraction(wn_random)
            stat[i] = kurtosis(wn_random) 
        }

        return(abs(median(stat) - kurtosis(wn_data)))
    }

    #' Fit alpha to a slip_mat, given it's corner wavenumber parameters
    #' We do this by optimizing 'alpha_goodness_of_fit'
    #' 
    #' @param slip_mat the slip matrix to fix
    #' @param wavenumber_par a length 2 vector with the wavenumber par for the slip
    #'        matrix. These need to be estimated, but our checks suggest they are often
    #'        insensitive to alpha
    #' @return The best alpha value
    #'
    fit_alpha<-function(slip_mat, wavenumber_par){

        fitted_alpha = optimize(alpha_goodness_of_fit, c(1.0e-06, 2.0), 
            slip_mat = slip_mat, wavenumber_par = wavenumber_par, 
            nreps=200)$minimum

        return(fitted_alpha)
    }


    #' Test case -- this methods seems ok, if not perfect, and it makes intuitive sense
    test_fit_alpha<-function(){

        test_alphas = c(1.1, 1.2, 1.5, 1.8, 1.9)

        for(alpha in test_alphas){

            print(paste0('Testing alpha = ', alpha))
            
            fitted_alphas=rep(NA, 100)

            # Simulate many SFFM, fit alpha to each
            for(i in 1:length(fitted_alphas)){

                #print(i)

                # Wavenumbers from raster[1]
                wavenumber_par = c(0.0623, 0.0771)

                sim_mat = simulate_sffm(c(wavenumber_par, alpha), as.matrix(rod$myrasts[[1]]))

                fitted_alphas[i] = fit_alpha(sim_mat, wavenumber_par)
            }

            # The median of these should be an ok estimator (if not highly accurate)
            print('Summary')
            print(summary(fitted_alphas))
            print('SD:')
            print(sd(fitted_alphas))

            median_CI = quantile_CI(fitted_alphas, level=0.99)
            print('Median 99% confidence interval')
            print(median_CI)

            if( (median_CI$ci[1]<= alpha) & (median_CI$ci[2]>=alpha) ){
                print('PASS')
            }else{
                print('FAIL')
            }
        }

    }  


    fit_all_FFI_alpha<-function(){

        fitted_alphas = rep(NA, length(rod$myrasts))

        for(i in 1:length(rod$myrasts)){

            print(i)

            ffi = as.matrix(rod$myrasts[[i]])

            # FFI parameters in numerical space
            ffi_par = c(rod$p1$kcxN[i], rod$p1$kcyN[i])

            fitted_alphas[i] = fit_alpha(ffi, ffi_par)

        }

        return(environment())
    }


    #' Compute a quantile on the upper 80% of the slip
    quantile_plus80<-function(slip_matrix, p=0.5){
        # Compute the median slip on the that part of the slip containing
        # 80% of the integrated slip
        sorted_slip = sort(c(slip_matrix))
        cumsum_sorted_slip = cumsum(sorted_slip)
        l = length(cumsum_sorted_slip)
        integrated_slip = cumsum_sorted_slip[l]
        keep = which(cumsum_sorted_slip > 0.2*integrated_slip)
        return(quantile(sorted_slip[keep], p=p))
    }


    ##' Code for a very large test of alpha fitting
    #massive_fitting_check<-function(parallel_fits){
    #    Nrast = 10
    #    best_alpha = rep(NA, Nrast)
    #    stat_store = vector("list", length=Nrast) #list()
    #    synthetic_mats_store = vector("list", length=Nrast) #list()
    #    for(n in 1:Nrast){
    #        print(paste0('n: ', n))
    #        # Generate 50 synthetic slips for each value of alpha
    #        synthetic_mats_store[[n]] = vector('list', length=11) #list()
    #        stat_store[[n]] = vector('list', length=11) #list()
    #        slip_mat = as.matrix(rod$myrasts[[n]])
    #        for(i in 1:11){
    #            print(i)
    #            synthetic_mats_store[[n]][[i]] = vector('list', length=50) #list()
    #            stat_store[[n]][[i]] = vector('list', length=50) #list()
    #            for(j in 1:50){
    #                synthetic_mats_store[[n]][[i]][[j]] = 
    #                    simulate_sffm(c(parallel_fits[[n]]$best_pars[[i]]$par, alphas[i]), slip_mat)
    #                stat_store[[n]][[i]][[j]] = 
    #                    fit_alpha(synthetic_mats_store[[n]][[i]][[j]], parallel_fits[[n]]$best_pars[[i]]$par)
    #            }
    #        }
    #    }

    #    return(environment())
    #}

    #X = best_alpha2(parallel_fits)

    #results_storage[[stable_dist_run]] = fit_all_FFI_alpha()


    # Test cases
    test_fit_alpha()

    # Fit FFI
    results_storage = fit_all_FFI_alpha()

    print('Summary of fits to all FFI')
    print(summary(results_storage$fitted_alphas))

    sink()

    return(environment())
}

alpha_results = list()
for(i in 1:4){
    # Reproducible
    set.seed(1234)

    alpha_results[[i]] = test_alpha_fit_and_fit_to_ffi(stable_dist_runs[i])
}

# Get summary of fitted alpha values for each FFI and each SFFM
#print(lapply(results_storage, f<-function(x) summary(x$fitted_alphas)))
