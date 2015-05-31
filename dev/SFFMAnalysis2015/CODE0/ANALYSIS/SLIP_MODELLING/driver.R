###############################################################################
# driver routine for analysis. 
#
# AUTHORS
# Gareth Davies, Geoscience Australia, 2013-2015,
# gareth.davies.ga.code@gmail.com
#
###############################################################################
library(methods) # For Rscript
library(parallel)

# Define models to run
configuration_data = list(

   S_GA_ = list(
                fourierFun = '2parSom', 
                noise_distribution = 'gaussian',
                negative_slip_removal_method = 'abs',
                spatial_slip_decay = 'none',
                RECENTRE_SLIP=TRUE,
                hurst=1.0,
                use_test_data=FALSE),
   
    S_GAF = list(
                fourierFun = '2parSom', 
                noise_distribution = 'gaussian',
                negative_slip_removal_method = 'abs',
                spatial_slip_decay = 'gaussian',
                RECENTRE_SLIP=TRUE,
                hurst=1.0,
                use_test_data=FALSE),
   
    S_GC_ = list(
                fourierFun = '2parSom', 
                noise_distribution = 'gaussian',
                negative_slip_removal_method = 'clip',
                spatial_slip_decay = 'none',
                RECENTRE_SLIP=TRUE,
                hurst=1.0,
                use_test_data=FALSE),
   
    S_GCF = list(
                fourierFun = '2parSom', 
                noise_distribution = 'gaussian',
                negative_slip_removal_method = 'clip',
                spatial_slip_decay = 'gaussian',
                RECENTRE_SLIP=TRUE,
                hurst=1.0,
                use_test_data=FALSE),

   S_SA_ = list(
                fourierFun = '2parSom', 
                noise_distribution = 'stable',
                negative_slip_removal_method = 'abs',
                spatial_slip_decay = 'none',
                RECENTRE_SLIP=TRUE,
                hurst=1.0,
                use_test_data=FALSE),
   
    S_SAF = list(
                fourierFun = '2parSom', 
                noise_distribution = 'stable',
                negative_slip_removal_method = 'abs',
                spatial_slip_decay = 'gaussian',
                RECENTRE_SLIP=TRUE,
                hurst=1.0,
                use_test_data=FALSE),
   
    S_SC_ = list(
                fourierFun = '2parSom', 
                noise_distribution = 'stable',
                negative_slip_removal_method = 'clip',
                spatial_slip_decay = 'none',
                RECENTRE_SLIP=TRUE,
                hurst=1.0,
                use_test_data=FALSE),
   
    S_SCF = list(
                fourierFun = '2parSom', 
                noise_distribution = 'stable',
                negative_slip_removal_method = 'clip',
                spatial_slip_decay = 'gaussian',
                RECENTRE_SLIP=TRUE,
                hurst=1.0,
                use_test_data=FALSE)
    
    )
                       



###############################################################################

# Ensure that names in the configuration_data are unique
configuration_names = names(configuration_data)
stopifnot(length(unique(configuration_names)) == length(configuration_data))

# Main loop -- repeat analysis for multiple model configuration setups
for(config_name in configuration_names){

    print(paste0('Running with configuration_data ', config_name))

    # Get global parameters + configuration_data as an environment
    # This is assumed to be in the search path of many routines
    config_pars = list2env(configuration_data[[config_name]])
    print(as.list(config_pars))

    print('Getting static config pars...')
    source('static_configuration_parameters.R', local = config_pars)
    print('....Done')

    # Give the 'same' random numbers to the master in each iteration of the loop
    # A side effect is that similar e.g. random phase parameters are generated
    # for different config_pars, which makes comparison easier
    initial_random_seed = .Random.seed
    

    # ID flag to add to file names etc
    runStamp = paste0(round(as.numeric(Sys.time())*100000), 
                      '_', config_name, 
                      '_NST_', config_pars$negative_slip_removal_method,
                      '_SSD_', config_pars$spatial_slip_decay,
                      '_RCS_', config_pars$RECENTRE_SLIP,
                      '_hurst_', config_pars$hurst,
                      '_noise_', config_pars$noise_distribution)

    # Ensure output directories exist
    dir.create(config_pars$quickfig_folder, showWarnings = FALSE, 
        recursive = TRUE)
    dir.create(config_pars$pdf_folder, showWarnings = FALSE, 
        recursive = TRUE)
    dir.create(config_pars$rimages_folder, showWarnings = FALSE, 
        recursive = TRUE)
    dir.create(config_pars$deformation_folder, showWarnings = FALSE, 
        recursive = TRUE)

    # Name for rimage file (repeatedly back-up the R session here)
    r_image_output_file = paste0(config_pars$rimages_folder, 
        '/run_results_', runStamp, '.Rdata')

    # Get the main analysis routines (which rely on config_pars)
    source('high_level_analysis.R')
    source('stable.R')

    print('Export everything to the cluster...')
    clusterExport(cl = config_pars$cluster, varlist = ls(), envir = environment())
    print('.....Done')

    # Fit the SRCMOD data to the model defined by config_pars
    print('Running run_on_data ...')
    rod = run_on_data(runStamp)

    # Save the session
    save.image(r_image_output_file)
    
    # For every fitted model, check that we can accurately fit a model like
    # that

    if(config_pars$num_test_fit_sims > 0){
        print('Simulation testing ...')
        allPars = lapply(rod$all_fits, myfun<-function(x) x$par)
        #testAllFit = mapply(test_statistical_fit2, allPars, rod$myrasts, 
        #    nsim=config_pars$num_test_fit_sims, SIMPLIFY=FALSE)

        testAllFit = clusterMap(
            cl = config_pars$cluster,
            test_statistical_fit2, 
            allPars, 
            rod$myrasts, 
            nsim=config_pars$num_test_fit_sims, 
            SIMPLIFY=FALSE,
            .scheduling='dynamic')

        test_fit_summary(testAllFit, rod$myrasts, runStamp)
    }

    # Plot the simulated raster kcx values for each finite fault
    plot_simulated_raster_kcx_kcy_Values(rod, runStamp)

    # Backup the session
    save.image(r_image_output_file)

    if(config_pars$num_synth_tsunami>0){
        # Simulate tsunami from synthetic FFMs with parameters selected from the
        # SRCMOD data
        print('Running simulate_many_tsunami ...')
        tsuMod = simulate_many_tsunami(rod, Nsim = config_pars$num_synth_tsunami,
            runStamp)

        # Backup the session
        save.image(r_image_output_file)

        #
        # Some quick plots for checking -- could do this in parallel
        #

        print('Plotting ...')
        # Plot synthetic deformation points
        synthetic_deformation_folders = 
            dir(paste0(config_pars$deformation_folder,'/',runStamp), 
                full.names = TRUE)
        clusterMap(cl=config_pars$cluster, plot_def_folder, 
            defFolder=as.list(synthetic_deformation_folders),
            .scheduling='dynamic')
        #for(i in synthetic_deformation_folders){
        #    print(i)
        #    plot_def_folder(i, cex = 2)
        #}

        # Eyeball Okada Deformation. All one file, should be serial
        pdf(paste0(config_pars$pdf_folder,'/OkadaDef_', runStamp, '.pdf'), 
            width = 6, height = 5)
        for(i in 1:length(tsuMod)){
            if(length(tsuMod[[i]]) == 1) next
            plot_okadaDef(tsuMod[[i]]$OkadaDef)
            title(paste0('STRK: ',  rod$myrastMetaData$STRK[i]), outer = TRUE, 
                line = -2)
            title(i, outer = TRUE, line = -3)
        }
        dev.off()
    }

    # Get 'Equivalent uniform slip earthquake parameters'
    write_equivalent_uniform_slip_eq_pars(rod, runStamp)

    # Clear the major memory hogs
    rm(list=c('rod', 'tsuMod', 'testAllFit'))
    # Stop the cluster (seem to get bugs if we don't do this)
    stopCluster(cl = config_pars$cluster)
    # Remove config_pars (we will read everything again in the loop)
    rm('config_pars')

}
