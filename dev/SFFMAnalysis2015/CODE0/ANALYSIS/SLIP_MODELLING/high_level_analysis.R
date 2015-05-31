###############################################################################
#
# High-level main analysis. 
#
# AUTHORS
# Gareth Davies, Geoscience Australia, 2013-2015, gareth.davies.ga.code@gmail.com
#
###############################################################################

library(parallel)
suppressPackageStartupMessages(library(raster))

# Key fitting routines
source('fit_simulate_earthquake.R', local=TRUE)

# We must have config_pars
if (!exists('config_pars')){
    stop('config_pars must exist for high_level_analysis.R to be used')
}

###############################################################################
#'
#' Major analysis routine, wrapped in a function which returns its environment
#'
#' @param runStamp is a character string which is included in filenames as an
#'        ID
#' @return the function environment
run_on_data<-function(runStamp){

    myrastlist = config_pars$myrastlist
    myrastMetaData = config_pars$myrastMetaData

    # Read rasters
    myrasts = list()
    myrast_stats = list()
    print('    Reading rasters ...')
    for(i in 1:length(myrastlist)){
        myrasts[[i]] = raster(myrastlist[i])
        myrast_stats[[i]] = Useful_eq_stats(list(myrasts[[i]]))
    }

    # Extract Mw from the filename
    myrasts_Mw = myrastMetaData$Mw


    # Use the cluster, checking that random seeds are unique
    config_pars$check_unique_random_seeds()

    # Fit the statistical model
    print('    Fitting rasters in parallel ...')
    #all_fits = mclapply(myrasts, fit_slip_parameters, mc.preschedule=FALSE, 
    #    mc.cores=config_pars$MC_CORES)
    all_fits = clusterMap(
        cl = config_pars$cluster, 
        fun = fit_slip_parameters, 
        m1 = myrasts, 
        .scheduling='dynamic')

    config_pars$check_unique_random_seeds()

    print('    .... Finished fit')

    p1 = get_fitpar(all_fits, myrasts)

    # Find regression relations between Mw and log10(kcx), Mw and log10(kcy)
    myrast_kcx = p1$kcx #p1$Kcx/myrastMetaData$LEN
    myrast_kcy = p1$kcy #p1$Kcy/myrastMetaData$WID
  
    if(config_pars$self_similar_kcx_mw_regression){
        # 'Self-similar' regression with coefficients of 0.5 
        lmx = lm(log10(myrast_kcx) ~ offset(-0.5*myrasts_Mw))
        lmy = lm(log10(myrast_kcy) ~ offset(-0.5*myrasts_Mw))
    }else{
        # Raw regression 
        lmx = lm(log10(myrast_kcx) ~ myrasts_Mw)
        lmy = lm(log10(myrast_kcy) ~ myrasts_Mw)
    }


    # Scatterplot of kcx-residual vs kcy-residual
    png(paste0(config_pars$quickfig_folder,'/kcx_kcy_regression_resid_', 
            runStamp, '.png'),
        width = 10, height = 6, units = 'in', res = 300)
    plot(lmx$resid, lmy$resid, 
        main='Scatterplot of Residuals after predicting log10(kcx) and log10(kcy) against Mw')
    dev.off()
   
    # Plot regression diagnostics kcx 
    png(paste0(config_pars$quickfig_folder,'/kcx_regression_diagnostics_', 
            runStamp, '.png'),
        width = 10, height = 6, units = 'in', res = 300)
    par(mfrow=c(2,2))
    plot(lmx)
    dev.off()
    
    # Plot regression diagnostics kcy
    png(paste0(config_pars$quickfig_folder,'/kcy_regression_diagnostics_', 
            runStamp, '.png'),
        width = 10, height = 6, units = 'in', res = 300)
    par(mfrow=c(2,2))
    plot(lmy)
    dev.off()

    # Generate synthetic rasters with simulation
    Nsim = config_pars$num_sims_earthquake_statistics

    if(RNGkind()[1] != "L'Ecuyer-CMRG"){
        stop("Must use RNG L'Ecuyer-CMRG for mcmapply with random numbers")
    }

    # Update the .Random.seed on the master (to ensure that if this function is
    # called in a serial loop, the seed changes)
    #.Random.seed = change_random_seed()

    #sim_rasts_and_par = mcmapply(
    #    synthetic_slip_from_fitted_slip, 
    #    all_fits, 
    #    myrasts, 
    #    myrasts_Mw, 
    #    Nsim, 
    #    list(lmx), 
    #    list(lmy),
    #    config_pars$random_kcx_kcy_simulations,
    #    SIMPLIFY=FALSE,
    #    mc.cores=config_pars$MC_CORES,
    #    mc.preschedule=FALSE)

    config_pars$check_unique_random_seeds()
    print('    Simulating rasters in parallel ...')

    sim_rasts_and_par = clusterMap(
        cl = config_pars$cluster,
        fun=synthetic_slip_from_fitted_slip,
        all_fits,
        myrasts,
        myrasts_Mw,
        Nsim,
        list(lmx),
        list(lmy),
        config_pars$random_kcx_kcy_simulations,
        .scheduling='dynamic')

    config_pars$check_unique_random_seeds()
    print('     ......Finished')
        
        
    # Extract list of simulated raster sets
    sim_rasts = lapply(sim_rasts_and_par, f<-function(x) x[[1]])

    # Extract parameters for simulated raster sets
    sim_rast_Pars = lapply(sim_rasts_and_par, f<-function(x) x[[2]])
    rm(sim_rasts_and_par)  # Free memory

    # Compute statistics from simulated rasters
    config_pars$check_unique_random_seeds()
    print('    Computing stats of sim_rasts ...')

    #myrast_simstats=mclapply(
    #    sim_rasts,
    #    Useful_eq_stats,
    #    mc.cores=config_pars$MC_CORES,
    #    mc.preschedule=FALSE)
    myrast_simstats = clusterMap(
        cl = config_pars$cluster,
        fun=Useful_eq_stats, 
        sim_rasts,
        .scheduling='dynamic')
    config_pars$check_unique_random_seeds()
    print('    ......Finished')

    print(    'Plotting...')

    # Boxplot the simstats and compare with the data
    Nevents = length(myrasts)
    mymat = matrix(NA, ncol=Nevents, nrow=Nsim)
    for(i in 1:Nevents){
         mymat[,i] = myrast_simstats[[i]]$median_80_slip/
            myrast_simstats[[i]]$quant_95_80_slip
     }

    png(paste0(config_pars$quickfig_folder, '/Sim_meas_skew_', runStamp,'.png'),
        width=10, height=6, units='in', res=300)

        boxplot(mymat, ylim=c(0.1,0.9), 
            main='Simulated and measured "skewness" stat')

        measured_stat = unlist(lapply(
            myrast_stats, 
            myfun<-function(x) x$median_80_slip/x$quant_95_80_slip ))
        points(1:Nevents, measured_stat, t='p', col=2, pch=19)

    dev.off()

    ###########################################################################
    # Make a multi-page plot to help visually compare the simulated + real
    # models for goodness-of-fit. 
    #
    transGrey = rgb(100, 100, 100, alpha=50, maxColorValue=255)
    pdf(paste0(config_pars$pdf_folder, '/SRCMOD_vs_simulated_', runStamp,'.pdf'), 
        width=15, height=9)
    par(mfrow=c(2,3))
    for(i in 1:Nevents){
        print(paste0('    ', i))
        ammi = as.matrix(myrasts[[i]]) # Just coerce to matrix once
        m1 = Mod(fft(ammi))
        tmp = get_wavenumbers(myrasts[[i]])

        # PLOT 1 -- sorted slip
        plot(sort(ammi), 
             main=paste('Sorted slip',i), xlab='Cell number', ylab='Slip',
                    ylim=c(0, max(ammi)))

        # Simulate models, and plot their sorted slip 
        Nreps = config_pars$num_sims_earthquake_statistics

        for(j in 1:Nreps){
            tmpfit = sim_rasts[[i]][[j]] 
            sorted_slip = sort(as.matrix(tmpfit))
            points(sorted_slip, t='l', col=transGrey)
        }
        legend('topleft', legend=c('Data', 'Model Realisations'), pch=c(1,NA), 
            lty=c(NA,1), col=c('black', 'grey'))

        # PLOT 2 -- spectrum
        plot(sqrt(tmp[[1]]**2 + tmp[[2]]**2), m1, 
            main=paste('Fourier Spectrum', i), xlab=expression(k[r]), 
                ylab='Spectrum')

        Fitted_wav = tmp[[1]]*0
        for(j in 1:Nreps){
            tmpfit = sim_rasts[[i]][[j]]
            amtmp = as.matrix(tmpfit)                      
            Fitted_wav = Fitted_wav+Mod(fft(amtmp))
            points(sqrt(tmp[[1]]**2 + tmp[[2]]**2), 
                Mod(fft(amtmp)), t='p',col=transGrey, pch='.')
        }

        Fitted_wav = Fitted_wav/Nreps
        points(sqrt(tmp[[1]]**2 + tmp[[2]]**2), Fitted_wav, t='p', col=2,
            pch=19, cex=0.5)
        legend('topright', legend=c('Data', 'Mean of Models'), 
            col=c('black', 'red'), pch=c(1,19))


        # PLOTS 3 4 5 6
        #
        # Make real/random slip surface plots
        # Ensure they have the same zlim
        randMod1 = sim_rasts[[i]][[1]] 
        randMod2 = sim_rasts[[i]][[2]] 
        randMod3 = sim_rasts[[i]][[3]] 
        ZLIM = c(0,max(max(ammi),
                       max(as.matrix(randMod1)), 
                       max(as.matrix(randMod2)), 
                       max(as.matrix(randMod3)) 
                      )
            )

        plot(myrasts[[i]], asp=1, zlim=ZLIM)
        title(basename(myrasts[[i]]@file@name))

        plot(randMod1, asp=1, zlim=ZLIM)
        title('A random realisation')
        
        plot(randMod2, asp=1, zlim=ZLIM)
        title('A random realisation')

        plot(randMod3, asp=1, zlim=ZLIM)
        title('A random realisation')

    }
    dev.off()

    ###########################################################################
    # OTHER PLOTS

    # Any temporal patterns in Kcx??
    pdf(paste0(config_pars$quickfig_folder, '/random_questions_', runStamp,
        '.pdf'))

        plot(p1$Kcx, p1$Kcy, main=paste('Is', expression(K[cx] == K[cy]), "?"))
        abline(0,1)

        barplot(p1$Kcx, main='Does Kcx increase over time??')

        # Any relation with Kcx + Mw??
        myMagCat = round(myrasts_Mw*2)/2
        myMagCatUnique = sort(unique(myMagCat))
        plot(p1$Kcx, p1$Kcy, col=myMagCat*2-12)
        legend('topleft', legend=as.character(myMagCatUnique), 
            col=myMagCatUnique*2-12, pch=1)

    dev.off()

    # Is Kcx in any sense 'stable' for a single event?
    eventCode=paste(myrastMetaData$Date, myrastMetaData$Location)
    png(paste0(config_pars$quickfig_folder, '/Kcx_events', runStamp,'.png'), 
            width=17, height=10, units='in', res=300)

        par(family='serif')
        par(mar=c(20,4,4,1))
        boxplot(p1$Kcx ~ as.factor(eventCode), las=2)
        title(expression(paste('Variation in ', K[cx], 
            ' for different FFMs of the same event')), cex.main=2)
        points(as.factor(eventCode),p1$Kcx)

    dev.off()

    # Variability in kcx for the same event 
    png(paste0(config_pars$quickfig_folder, '/kcx_events', runStamp, '.png'),
        width=17, height=10, 
        units='in', res=300)

        myrast_kcx = p1$Kcx/myrastMetaData$LEN
        par(mar=c(20,4,4,1))
        boxplot(log10(myrast_kcx) ~ as.factor(eventCode), las=2)
        title(expression(paste('Variation in ', log[10](k[cx]), 
            ' for different FFMs of the same event')), cex.main=2)
        points(as.factor(eventCode), log10(myrast_kcx))

    dev.off()

    # Try to distinguish single events
    png(paste0(config_pars$quickfig_folder, '/Mw_kcx', runStamp, '.png'),
        width=10, height=10, res=300, units='in')

        plotCodes = c(letters, LETTERS)
        plot(myrasts_Mw, log10(myrast_kcx), 
            pch=plotCodes[as.numeric(as.factor(eventCode))])

    dev.off()

    # Look at the 'high slip' in the synthetic models and the FFMs
    # NOTE: This does not account for the fraction of the rupture containing
    # most seismic moment, like we did in the paper
    max_95th = unlist(lapply(myrast_simstats, 
        f<-function(x) max(x$quant95_slip)))
    min_95th = unlist(lapply(myrast_simstats, 
        f<-function(x) min(x$quant95_slip)))
    median_95th = unlist(lapply(myrast_simstats,
        f<-function(x) median(x$quant95_slip)))
    ffm_95th = unlist(lapply(myrasts, 
        f<-function(x) quantile(as.matrix(x),0.95)))

    png(paste0(config_pars$quickfig_folder, '/HighSlip_biases', runStamp,'.png'),
        width=10, height=10, 
        res=300, units='in')

        YLIM = c(min(c(min_95th, ffm_95th)), max(c(max_95th, ffm_95th)))
        plot(ffm_95th, median_95th, log='xy', 
            main='95th percentile of slip in FFMs vs synthetic models',
            pch=round(myrasts_Mw), ylim=YLIM)
        arrows(ffm_95th, min_95th, ffm_95th, max_95th, length=0.)
        points(ffm_95th, max_95th, col=2, cex=1.0, pch=19)
        points(ffm_95th, min_95th, col=3, cex=1.0, pch=19)
        grid()
        abline(0, 1)

    dev.off()

    max_50th = unlist(lapply(myrast_simstats, 
        f<-function(x) max(x$median_slip)))
    min_50th = unlist(lapply(myrast_simstats, 
        f<-function(x) min(x$median_slip)))
    median_50th = unlist(lapply(myrast_simstats, 
        f<-function(x) median(x$median_slip)))
    ffm_50th = unlist(lapply(myrasts,
        f<-function(x) quantile(as.matrix(x),0.5)))

    png(paste0(config_pars$quickfig_folder, '/MedianSlip_biases', runStamp,'.png'), width=10, 
        height=10, res=300, units='in')

        YLIM = c(min(c(min_50th, ffm_50th)), max(c(max_50th, ffm_50th)))
        plot(ffm_50th, median_50th, log='xy',
            main='Median slip in FFMs vs synthetic models', 
            pch=round(myrasts_Mw))
        arrows(ffm_50th, min_50th, ffm_50th, max_50th, length=0.)
        points(ffm_50th, max_50th, col=2, cex=1.0, pch=19)
        points(ffm_50th, min_50th, col=3, cex=1.0, pch=19)
        grid()
        abline(0, 1)

    dev.off()

    return(environment())
}   

###############################################################################
#'
#' Given some best fit corner wavenumbers and a raster, generate a number of
#' synthetic rasters with the same parameters
#'
#' @param modelFit result of call to fit_slip_parameters, containing model fits
#' @param rast template raster for simulation (used to determine dx, dy, etc)
#' @param Mw rupture magnitude, used to predict parameters
#' @param Nsim Number of synthetic slip surfaces to make
#' @param lmx Linear model output for the kcx parameter
#' @param lmy Linear model output for the kcy parameter
#' @param random_kcx_kcy_simulations TRUE/FALSE If true use the magnitude +
#'        linear models to predict new kcx/kcy including a random component
#'        from the regression fit. Otherwise, just use modelFit$par 
#' @return list of length 2, containing a list with all the synthetic rasters,
#'          and a matrix of their synthetic parameters
#'
synthetic_slip_from_fitted_slip<- function(
    modelFit,
    rast,
    Mw,
    Nsim,
    lmx,
    lmy,
    random_kcx_kcy_simulations){

    # Storage lists
    sim_rasts = as.list(rep(NA,Nsim))
    simParStore = matrix(NA,ncol=length(modelFit$par),nrow=Nsim)

    for(j in 1:Nsim){

        # Optionally use simulated kcx/kcy from regression
        if(!random_kcx_kcy_simulations){
           simPar=modelFit$par

        }else{
            # Assume first model parameter is kcx, second is kcy
            simPar=modelFit$par*NA

            # Simulate new kcx/kcy from regression relations
         
            ## APPROACH 1: 
            # logMean=predict(lmx,data=data.frame(myrasts_Mw=Mw))
            # logSD=summary(lmx)$sigma
            ## Compute random kcx, and convert to numerical space by
            ## multiplication with raster res 
            # simPar[1] = (10**(rnorm(1,logMean,logSD)))*res(rast)[1]

            # logMean=predict(lmy,data=data.frame(myrasts_Mw=Mw))
            # logSD=summary(lmy)$sigma
            ## Compute random kcx, and convert to numerical space by
            ## multiplication with raster res 
            # simPar[2] = (10**(rnorm(1,logMean,logSD)))*res(rast)[2]

            # APPROACH 2: RESAMPLE RESIDUALS FROM REGRESSION
            # Still reasonable with correlation in residuals of kcx,kcy when
            # predicted using Mw. Use bootstrap sampling of residuals

            myint = sample(1:length(lmx$resid), size=1)

            lmx_prediction = predict(lmx, 
                newdata=data.frame(myrasts_Mw=Mw)) + 
                lmx$resid[myint]
            simPar[1] = (10**(lmx_prediction))*res(rast)[1]

            lmy_prediction = predict(lmy, 
                newdata=data.frame(myrasts_Mw=Mw)) +
                lmy$resid[myint]
            simPar[2]=(10**(lmy_prediction))*res(rast)[2]

        }

        simParStore[j,] = simPar
        sim_rasts[[j]] = simulate_sffm(simPar,rast)
    }
    return(list(sim_rasts, simParStore))
}

###############################################################################
#'
#' Utility to extract parameters from a list of model fits and a list of their
#' rasters
#'
#' @param all_fits list of objects from optim containing model fits
#' @param myrasts list of rasters corresponding to all_fits
#' @return the function environment, which contains various useful parameters
#' 
get_fitpar<-function(all_fits, myrasts){

    # n parameter
    nn=unlist(lapply(all_fits, myfun<-function(x) x$par[3]))

    # Numerical kcx,kcy
    kcxN=unlist(lapply(all_fits, myfun<-function(x) x$par[1]))
    kcyN=unlist(lapply(all_fits, myfun<-function(x) x$par[2]))

    # Number of cells in raster x,y
    rastNX=unlist(lapply(myrasts, myfun<-function(x) dim(x)[2]))
    rastNY=unlist(lapply(myrasts, myfun<-function(x) dim(x)[1]))

    #browser()

    stopifnot(length(nn) == length(kcxN))
    stopifnot(length(kcxN) == length(kcyN))
    stopifnot(length(kcyN) == length(rastNX))
    stopifnot(length(kcyN) == length(rastNY))

    # Kcx parameter
    Kcx=kcxN*rastNX
    Kcy=kcyN*rastNY

    # Convergence flag
    converged=unlist(lapply(all_fits, myfun<-function(x) x$convergence))

    # Physical kcx,kcy
    kcx=kcxN/unlist(lapply(myrasts,myfun<-function(x) res(x)[1]))
    kcy=kcyN/unlist(lapply(myrasts,myfun<-function(x) res(x)[2]))
    
    return(environment())
}


###############################################################################
#'
#' Code to simulate 'random' tsunami deformation, based on FFM
#'
#' @param myrast_FSP name of FSP file corresponding to the FFM (it contains
#'        info not in the raster itself)
#' @param myrast raster corresponding to the FFM
#' @param reg_pars fitted model parameters
#' @param strike rupture strike 
#' @param dip rupture dip
#' @param rake rupture rake (NOTE: only used if Rake is not a column in the FSP
#'             file)
#' @param nsim number of synthetic events to simulate
#' @param runStamp runStamp ID for the simulation
#' @param simSlipPars parameters for the synthetic events. If NULL, use reg_pars
#' 
#' @return a list with lots of important info
#'
simulate_tsunami<-function(
    myrast_FSP, 
    myrast, 
    reg_pars, 
    strike, 
    dip,
    rake, 
    nsim, 
    runStamp, 
    simSlipPars){

    
    #test=simulate_tsunami(myrast_FSP=myrastFSPnames[55], myrast=rod$myrasts[[55]], 
    #    reg_pars=rod$all_fits[[55]]$par,strike=rod$myrastMetaData$STRK[55],
    #    dip=rod$myrastMetaData$DIP[55], rake=rod$myrastMetaData$RAKE[55], return_xyDef=TRUE, 
    #    nsim=2,
    #    runStamp=runStamp, simSlipPars=rod$sim_Rast_Pars[[55]])
    

    # Run the 'real' FFM deformation
    myDeformation = try(
        runOkada(
            myrast_FSP, 
            myrast,
            strk=strike,
            dip=dip,
            rake=rake,
            return_xyDef=TRUE,
            simID="Orig", 
            runStamp=runStamp)
        )

    # If the 'try' catches an error, pass it to a list so the cluster
    # still handles it ok
    if(class(myDeformation)=='try-error'){
        return(list(myDeformation))
    }

    simSlip=as.list(rep(NA,nsim))
    simDef=as.list(rep(NA,nsim))

    # Make synthetic deformation, and run that
    if(nsim > 0){
        for(j in 1:nsim){

            # If simSlipPars is not NULL, it contains kx/ky parameters for the
            # simulated earthquakes. Otherwise, reg_pars from the function call
            # is used
            if(!is.null(simSlipPars)){
                reg_pars = simSlipPars[j,]
            } 

            simSlip[[j]] = simulate_sffm(reg_pars, myrast)

            simDef[[j]] = runOkada(
                myrast_FSP, 
                simSlip[[j]], 
                identical_rasters=FALSE,
                strk=strike,
                dip=dip,
                rake=rake,
                return_xyDef=TRUE,
                simID=as.character(j),
                runStamp=runStamp)
        }
    }

    return(list(fspFile=myrast_FSP, 
                OkadaDef=myDeformation,
                simSlip=simSlip,
                simDef=simDef))
}

###############################################################################
#'
#' Apply simulate_tsunami to rasters used in run_on_data, in parallel
#'
#' @param rod is the environment returned from run_on_data
#' @param Nsim value of 'nsim' passed to simulate_tsunami
#' @param runStamp runStamp ID for this run 
#'
#' @return the result of an mcmapply call to simulate_tsunami, applied over all
#'         the TIF files in rod$myrastlist
#'
simulate_many_tsunami<-function(rod, Nsim=1, runStamp='xxxxxxxxx', runparallel=TRUE){

    # FSP files for each raster
    myrast_FSPs = gsub('.tif', '.FSP', rod$myrastlist)
    
    # Extract regression paramaeter into list for mcmapply
    myRegPar = lapply(rod$all_fits, myfun<-function(x) x$par)

    # If running random kcx/kcy simulations, use the parameter values simulated
    # previously in rod [run_on_data() output]
    if(config_pars$random_kcx_kcy_simulations){
        mySimPar = rod$sim_rast_Pars
        if(rod$Nsim < Nsim){
            stop('ERROR: not enough random_kcx_kcy_simulations in rod to simulate this many tsunami')
        }
    }else{
        mySimPar=NULL
    }

    if(runparallel){ 
        library(parallel)

        if(RNGkind()[1] != "L'Ecuyer-CMRG"){
            stop("Must use RNG L'Ecuyer-CMRG for mcmapply with random numbers")
        }

        # Update the .Random.seed on the master (to ensure that if this function is
        # called in a loop, the seed changes)
        #.Random.seed = change_random_seed()

        #lotsOfTsunami=mcmapply(simulate_tsunami,
        #                       myrast_FSP=myrast_FSPs,
        #                       myrast=rod$myrasts, 
        #                       reg_pars=myRegPar,
        #                       strike=rod$myrastMetaData$STRK,
        #                       dip=rod$myrastMetaData$DIP,
        #                       rake=rod$myrastMetaData$RAKE,
        #                       nsim=rep(Nsim,length(myrast_FSPs)),
        #                       runStamp=runStamp,
        #                       simSlipPars=mySimPar, 
        #                       mc.preschedule=FALSE,
        #                       mc.cores=config_pars$MC_CORES,
        #                       SIMPLIFY=FALSE)
        config_pars$check_unique_random_seeds()

        lotsOfTsunami = clusterMap(
            cl = config_pars$cluster,
            fun=simulate_tsunami,
            myrast_FSP=myrast_FSPs,
            myrast=rod$myrasts, 
            reg_pars=myRegPar,
            strike=rod$myrastMetaData$STRK,
            dip=rod$myrastMetaData$DIP,
            rake=rod$myrastMetaData$RAKE,
            nsim=rep(Nsim,length(myrast_FSPs)),
            runStamp=runStamp,
            simSlipPars=mySimPar, 
            .scheduling='dynamic')

        config_pars$check_unique_random_seeds()

    }else{
        lotsOfTsunami = mapply(
            fun=simulate_tsunami,
            myrast_FSP=myrast_FSPs,
            myrast=rod$myrasts, 
            reg_pars=myRegPar,
            strike=rod$myrastMetaData$STRK,
            dip=rod$myrastMetaData$DIP,
            rake=rod$myrastMetaData$RAKE,
            nsim=rep(Nsim,length(myrast_FSPs)),
            runStamp=runStamp,
            simSlipPars=mySimPar)
    }
                           
    return(lotsOfTsunami)                          
}

###############################################################################
#test_statistical_fit<-function(
#    dataSource='../TEST_DATA/random_R/random_clipped_1_16x8/*.tif',
#    n=200){
#
#    # Now try on test rasters with Kcx=Kcy=1, and no spatial filtering / edge
#    # tapering etc. This should show the method works
#    myrast_testlist = Sys.glob(dataSource)[1:n]
#    myrasts_test = list()
#    for(i in 1:length(myrast_testlist)){
#        # NOTE: If the model includes distance decay, then the slip must be
#        # recentred. The simulated data was not recentred, so we force that here
#        myrasts_test[[i]]=recentre_slip(raster(myrast_testlist[i]))
#    }
#
#
#    if(RNGkind()[1] != "L'Ecuyer-CMRG"){
#        stop("Must use RNG L'Ecuyer-CMRG for mclapply with random numbers")
#    }
#
#    all_fits_test=mclapply(myrasts_test,fit_slip_parameters,mc.preschedule=FALSE,mc.cores=config_pars$MC_CORES)
#    p3=get_fitpar(all_fits_test,myrasts_test)
#
#    return(environment())
#}

#'
#' Simulate a model from reg_par, then apply fit_slip_parameters to the
#' simulated model
#'
#' Applied many times, this allows us to check the reliability of the
#' statistical techniques
#' 
#' @param reg_par model parameters (passed to simulate_sffm)
#' @param tg template slip model (passed to simulate_sffm)
#' @return result of fit_slip_parameters applied to the synthetic model
#'
randFitfun<-function(reg_par,tg=NULL){

    # Test whether we can back-fit reg_par on template raster tg
    suppressPackageStartupMessages(library(raster))

    m1 = simulate_sffm(reg_par, tg)

    output = fit_slip_parameters(m1)

    # Store the initial seed from the simulated model as an attribute (for
    # reproducibility)
    attr(output, 'initial_seed') = attr(m1, 'initial_seed')

    return(output)

}    

test_statistical_fit2<-function(reg_par, tg, nsim=100, runparallel=FALSE){
    # This can be used to test our ability to fit ANY model we can simulate,
    # with regression parameter vector reg_par, and template raster tg


    if(runparallel){
        if(RNGkind()[1] != "L'Ecuyer-CMRG"){
            stop("Must use RNG L'Ecuyer-CMRG for mclapply with random numbers")
        }

        # Update the .Random.seed on the master (to ensure that if this function is
        # called in a serial loop, the seed changes)
        #.Random.seed = change_random_seed()

        #output = mclapply( 
        #    rep(list(reg_par),nsim), 
        #    randFitfun,tg=tg, 
        #    mc.cores=config_pars$MC_CORES,
        #    mc.preschedule=FALSE)
        config_pars$check_unique_random_seeds()

        output = clusterMap(cl = config_pars$cluster,
            fun = randFitfun,
            rep(list(reg_par),nsim), 
            list(tg),
            .scheduling='dynamic')

        config_pars$check_unique_random_seeds()

    }else{

        output = lapply( 
            rep(list(reg_par),nsim), 
            randFitfun,
            tg=tg)

    }
 
    return(output)
}

###############################################################################
#' Extract info from a list of outputs from test_statistical_fit_2, which comes
#' from applying test_statistical_fit2 to many simulated rasters
#'
#' @param testAllFits list of outputs from test_statistical_fit_2
#' @param myrasts list of rasters corresponding to testAllFits
#' @param runStamp runstamp of this simulation
#' @return Nothing, but make various plots which show how the fitting performs
#'
test_fit_summary<-function(testAllFits,myrasts, runStamp){
    pdf(paste0(config_pars$pdf_folder, '/FitTest_', runStamp,'.pdf'))
    kcxFitMean=rep(NA,length(testAllFit))
    kcxFitSd=rep(NA,length(testAllFit))
    kcxTrue=rep(NA,length(testAllFit))
    kcyFitMean=rep(NA,length(testAllFit))
    kcyFitSd=rep(NA,length(testAllFit))
    kcyTrue=rep(NA,length(testAllFit))

    kcxAll=matrix(NA,ncol=length(testAllFits),nrow=length(testAllFits[[1]]))
    kcyAll=matrix(NA,ncol=length(testAllFits),nrow=length(testAllFits[[1]]))
    for(j in 1:length(testAllFit)){
        mm=matrix(unlist(lapply(testAllFit[[j]], f<-function(x) x$par)),ncol=2,byrow=T)
        par(mfrow=c(2,1))
        boxplot(mm)
        kcxFitMean[j]=mean(mm[,1])/res(myrasts[[j]])[1]
        kcxFitSd[j]=sd(mm[,1])/res(myrasts[[j]])[1]
        kcyFitMean[j]=mean(mm[,2])/res(myrasts[[j]])[2]
        kcyFitSd[j]=sd(mm[,2])/res(myrasts[[j]])[2]
        abline(h=rod$all_fits[[j]]$par/res(myrasts[[j]]),col=c('red','blue'))
        kcxTrue[j]=rod$all_fits[[j]]$par[1]/res(myrasts[[j]])[1]
        kcyTrue[j]=rod$all_fits[[j]]$par[2]/res(myrasts[[j]])[2]
        plot(rod$myrasts[[j]])

        kcxAll[,j]=mm[,1]/res(myrasts[[j]])[1]
        kcyAll[,j]=mm[,2]/res(myrasts[[j]])[2]
    }
    dev.off()

    # Plot the fitted vs true 
    png(paste0(config_pars$quickfig_folder, '/Fit_bias_', runStamp,'.png'), 
        width=12,height=8,units='in',res=300)
    par(family='serif')
    par(mfrow=c(2,2))
    plot(kcxFitMean,kcxTrue,main=expression(paste(k[cx],' fitted vs true')))
    abline(0,1); grid()
    plot(kcyFitMean,kcyTrue,main=expression(paste(k[cy],' fitted vs true')))
    abline(0,1); grid()
    plot(kcxFitMean, (kcxFitMean-kcxTrue)/kcxTrue, 
        main=expression(paste(k[cx],' relative error')))
    grid()
    plot(kcyFitMean, (kcyFitMean-kcyTrue)/kcyTrue, 
        main=expression(paste(k[cy], ' relative error')))
    grid()
    dev.off()

    #plot(kcxFitMean-kcxTrue, kcyFitMean-kcyTrue)
    png(paste0(config_pars$quickfig_folder, '/Fit_error_correlation_1_',
            runStamp,'.png'), 
        width=6,height=6,units='in',res=300)
    plot((kcxFitMean-kcxTrue)/kcxFitSd, (kcyFitMean-kcyTrue)/kcyFitSd)
    dev.off()

    png(paste0(config_pars$quickfig_folder, '/Fit_vs_true_',runStamp,'.png'), 
        width=12,height=9,units='in',res=300)
    par(mfrow=c(2,1))
    boxplot(kcxAll,main=expression(paste('Fitted ',k[cx])),log='y')
    points(kcxTrue,col='red',pch=19,cex=0.5); grid()
    boxplot(kcyAll,main=expression(paste('Fitted ',k[cy])),log='y')
    points(kcyTrue,col='red',pch=19,cex=0.5); grid()
    dev.off()
    return(NULL)
}

plot_tsunami_deformation_fourierspec<-function(tsuMod, runStamp){
    # Make plots of fourier spectrum vs radial wavenumber for okada deformation

    plotDIR=paste0(config_pars$quickfig_folder, '/okada_spectral_',runStamp)
    dir.create(plotDIR,showWarnings=FALSE,recursive=T)
    for(j in 1:length(tsuMod)){
        plotname=paste0(plotDIR,'/',
                       gsub('.FSP','.png',basename(names(tsuMod)[j]))
                )
        png(plotname,width=12,height=6,units='in',res=300)
        if(length(tsuMod[[j]])!=4){
            plot(c(0,1),c(0,1),main='FAILED TO PARSE')
            dev.off()
            next
        }
        # Get fourier spec / numerical wavenumbers
        m1=Mod(fft(as.matrix(tsuMod[[j]]$OkadaDef)))
        tmp=get_wavenumbers(m1)
        kx=tmp[[1]]
        ky=tmp[[2]]
        kr=sqrt(kx**2+ky**2)

        plot(kr,m1,log='y',pch='.',col='red')
        for(i in 1:50){
            m2=Mod(fft(as.matrix(tsuMod[[j]]$simDef[[i]])))
            points(kr,m2,pch='.',col=rgb(50,50,50,alpha=10,maxColorValue=255))
        }
        dev.off()
    }

}


compute_uniform_slipPars_from_FFMs<-function(Mw, dip, rake, peakSlipDepth, Slip_x_Area){
    # We are interested to compare tsunami from uniform-slip runs to tsunami
    # from FFMs and synthetic FFMs
    # To do this, we need an 'equivalent' uniform slip tsunami
    #
    # Here, the approach tries to be analgous to what a tsunami hazard study would do
    #
    # It has US_Mw = Mw
    #        US_dip = dip
    #        US_rake = rake
    #
    # Slip comes so that slip x Area in the FFM = slip x Area in the uniform slip model
    #
    # L, W come from the Blazer et al reverse fault relation
    #
    # The centroid location is determined by:
    #    1) Trying to make the FFM peak slip depth as the centroid depth
    #    2) If the fault would extend out of the earth, deepen it until it is 100m in the earth
    #

    # Use Blazer et al regression relations for length, width
    # NOTE: We use their ordinary least squares relations for reverse (thrust) faults, 
    #       To my knowledge OLS is the right approach for predicting log10(L)
    #       from Mw, even if there is error in Mw as well as log10(L). This
    #       is different to what they claim in the paper
    US_L=10**( -2.28+0.55*Mw) # Uniform Slip Length (km)
    US_W=10**( -1.86+0.46*Mw) # Uniform Slip Width (km)

    US_dip=dip
    US_rake=rake

    # Seismic Moment M0, in units dyne cm
    #M0=10**( (Mw+10.7)*(3./2.))
    #mu = 3.0e+11 # shear modulus, in dyne/cm^2

    A = US_L*US_W # Area in units km^2
    US_slip = Slip_x_Area/A # Slip, in units of m [note slip x Area is in units km^2]

    # Centroid Depth for uniform slip 
    US_centroidDepth=pmax(peakSlipDepth, 0.5*US_W*sin(dip/180*pi)+0.1)
    

    return(data.frame(Mw=Mw, length=US_L, width=US_W, 
                      dip=US_dip,rake=US_rake, slip=US_slip,
                      centroidDepth=US_centroidDepth))
}

Mw_2LW<-function(Mw) 10**( -c(2.28, 1.86)+c(0.55,0.46)*Mw)

##

write_equivalent_uniform_slip_eq_pars<-function(rod, runStamp){
    # Write 'equivalent' uniform slip parameters to a file
    # 
    # put here to keep 'runThis' clean

    # Mw
    FFM_Mw = rod$myrasts_Mw
    # Dip, rake in degrees
    FFM_dip = rod$myrastMetaData$DIP
    FFM_rake = rod$myrastMetaData$RAKE
    # down dip max depth, in km
    FFM_downDip_max = unlist(lapply(rod$myrast_stats, 
        f<-function(x) x$downDip_max_store))
    # Depth = down-dip-max*sin(dip) + depth of top of fault
    FFM_peakSlipDepth = FFM_downDip_max*sin(FFM_dip/180*pi) + 
        rod$myrastMetaData$Htop

    # slip x area. 
    # Units are m x (km^2)
    FFM_Slip_x_Area = rep(NA,length(rod$myrasts))
    for(i in 1:length(rod$myrasts)){
        FFM_Slip_x_Area[i] = sum(as.matrix(rod$myrasts[[i]]))*
           prod(res(rod$myrasts[[i]]))
    }

    uniformSlipPar = compute_uniform_slipPars_from_FFMs(FFM_Mw, FFM_dip, 
        FFM_rake, FFM_peakSlipDepth, FFM_Slip_x_Area)

    # Write to file, with enough meta-data to be sure of which event is which
    eqNumber = unlist(
        lapply(
            strsplit(gsub('.tif', '', basename(rod$myrastlist)), '_'), 
            f<-function(x) as.numeric(x[5]))
        )
    eqName = gsub('.tif', '', basename(rod$myrastlist))
    uniformSlipPar = cbind(data.frame(eqNumber=eqNumber, eqName=eqName), 
        uniformSlipPar)
    output_file_name = paste0(config_pars$deformation_folder, '/', runStamp, 
        '/FFM_uniformSlipPars.csv')
    dir.create(dirname(output_file_name), showWarnings = FALSE, 
        recursive = TRUE)
    write.table(uniformSlipPar, 
                file=output_file_name,
                sep=",", row.names=FALSE)

}

###

plot_simulated_raster_kcx_kcy_Values<-function(rod, runStamp){
    #@ Make plots of Mw vs boxplots with simulated kcx/kcy values, fitted
    #@ kcx/kcy values, and regression lines
    #
    #@ Put here to keep 'runThis' looking neat
    #
    #@ Note: rod$sim_rast_Pars contains NUMERICAL kcx,kcy, so I convert to physical in 'mapply'

    png(paste0(config_pars$quickfig_folder, '/Simulated_kcx_vs_Mw_', runStamp, '.png'), 
        width=8,height=6,units='in',res=300)
    par(family='serif')
    sim_kcx=mapply(f<-function(x,y) log10(x[,1]/res(y)[1]), rod$sim_rast_Pars, 
        rod$myrasts, SIMPLIFY=FALSE)
    boxplot(sim_kcx, at=rod$myrasts_Mw, boxwex=0.1,xaxt='n',ylim=c(-4,-1))
    title('Simulated FFM log10(kcx) values, vs Mw, with fitted values in red',
          xlab='Mw', ylab='log10(kcx)')
    axis(side=1)
    points(rod$myrasts_Mw,log10(rod$p1$kcx),col=2,pch=19)
    abline(coef(rod$lmx)[1], coef(rod$lmx)[2],col='blue')
    dev.off()

    png(paste0(config_pars$quickfig_folder, '/Simulated_kcy_vs_Mw_', runStamp, '.png'), 
        width=8,height=6,units='in',res=300)
    par(family='serif')
    sim_kcy=mapply(f<-function(x,y) log10(x[,2]/res(y)[2]), rod$sim_rast_Pars, 
        rod$myrasts, SIMPLIFY=FALSE)
    boxplot(sim_kcy, at=rod$myrasts_Mw, boxwex=0.1,xaxt='n',ylim=c(-4,-1))
    title('Simulated FFM log10(kcy) values, vs Mw, with fitted values in red',
          xlab='Mw', ylab='log10(kcy)')
    axis(side=1)
    points(rod$myrasts_Mw,log10(rod$p1$kcy),col=2,pch=19)
    abline(coef(rod$lmy)[1], coef(rod$lmy)[2],col='blue')
    dev.off()
}


## 

back_calculate_Mw<-function(myrasts){
    # Back-calculate Mw for all FFMs, assuming
    #  M0=10**( (Mw+10.7)*(3./2.))
    #  mu = 3.0e+11 # shear modulus, in dyne/cm^2
    lmr=length(myrasts)
    Mw_est=rep(NA,len=lmr)

    mu=3.0e+11

    for(i in 1:lmr){
        # Seismic moment, dyn/cm^2*cm*cm*cm = dyn*cm
        M0=mu*sum(as.matrix(myrasts[[i]]))*100.*prod(res(myrasts[[i]]))*(1000*100)*(1000*100)
        Mw_est[i] = log10(M0)*2/3 -10.7
    }
    return(Mw_est)
}

#' Quick check of how mean slip varies along-strike / down dip
plot_mean_slip_by_xy<-function(myrasts){
    
    pdf(paste0(config_pars$pdf_folder, '/down_dip_mean_slip.pdf'))

    for(i in 1:length(myrasts)){
        dm = slip_vs_distance(myrasts[[i]])
        plot(rowMeans(dm[[2]]), t='o')
        title(basename(myrasts[[i]]@file@name))
    }
    dev.off()

    pdf(paste0(config_pars$pdf_folder, '/along_strike_mean_slip.pdf'))
    for(i in 1:length(myrasts)){
        dm = slip_vs_distance(myrasts[[i]])
        plot(colMeans(dm[[2]]), t='o')
        title(basename(myrasts[[i]]@file@name))
    }
    dev.off()
}

#' Plots of the white noise to help identify likely good distributions
plot_white_noise<-function(myrasts, fit_dist='cauchy'){

    if(fit_dist == 'gumbel'){
        library(ismev)
    }else{
        # Assume fitdistrplus can handle it
        library(fitdistrplus)
        library(HyperbolicDist) # skew-laplace
        source('laplace.R')
    }

    store_fits = list()

    pdf(paste0(config_pars$pdf_folder, '/white_noise_', fit_dist,'.pdf'), width=10,height=8)

    for(i in 1:length(myrasts)){
        white_noise = white_noise_extraction(myrasts[[i]])
        par(mfrow=c(2,2))
        image(white_noise)
        title(paste0('white noise for \n', basename(myrasts[[i]]@file@name)))
        hist(white_noise)
        qqnorm(white_noise)
        qqline(white_noise)

        switch(fit_dist,
            'gumbel' = {
                # Fit a distribution
                gumbel_fit = gum.fit(c(white_noise))
                # Diagnostics
                plot(c(0,1), c(0,1), col=0, xlab="", ylab="")
                text(0.2, 0.7, paste0(signif(gumbel_fit$mle[1], 2), ', ', signif(gumbel_fit$mle[2], 2)))
                text(0.2, 0.3, paste0(signif(gumbel_fit$se[1], 2), ', ', signif(gumbel_fit$se[2], 2)))
                # Some plots
                par(mfrow=c(2,2))
                gum.diag(gumbel_fit)
            
                store_fits[[basename(myrasts[[i]]@file@name)]] = gumbel_fit
            },
                { # Assume fitdist can handle it

                if(fit_dist == 'lap'){
                    cauchy_fit = fitdist( c(white_noise), fit_dist, start = list(scale = mad(c(white_noise))) )
                    # Add some 'fake' parameters so the following plot codes work
                    #cauchy_fit$estimate = c(cauchy_fit$estimate, 0)
                    #cauchy_fit$sd = c(cauchy_fit$sd, 0)
                }else{
                    cauchy_fit = fitdist(c(white_noise), fit_dist)
                }
                #plot(c(0,1), c(0,1), col=0, xlab="", ylab="")
                #text(0.2, 0.7, paste0(signif(cauchy_fit$estimate[1], 2), ', ', signif(cauchy_fit$estimate[2], 2)))
                #text(0.2, 0.3, paste0(signif(cauchy_fit$sd[1], 2), ', ', signif(cauchy_fit$sd[2], 2)))
                # Some plots
                par(mfrow=c(2,2))
                #gum.diag(gumbel_fit)
                ppcomp(cauchy_fit)
                qqcomp(cauchy_fit)
                cdfcomp(cauchy_fit)
                denscomp(cauchy_fit)
                            
                store_fits[[basename(myrasts[[i]]@file@name)]] = cauchy_fit

            }            
        )
        
    }
    dev.off()

    return(store_fits)
}

#

#'
#' Quick plot of the spectrum vs a 1d wavenumber
plot_1d_spec<-function(myrast, fitted_par, add=FALSE, ...){
    wavenumber_matrices = get_wavenumbers(myrast)
    
    K_1d = ( (wavenumber_matrices[[1]]/fitted_par[1])**2 + 
             (wavenumber_matrices[[2]]/fitted_par[2])**2 )**0.5

    Fs = Mod(fft(as.matrix(myrast)))

    if(!add){
        plot(K_1d, Fs, log='y', ...)
    }else{
        points(K_1d, Fs, ...)
    }


    #Fs = c(Fs)
    #K_1d = c(K_1d)

    #fit = lm(log(Fs) ~ K_1d)
    #abline(coef(fit)[1], coef(fit)[2], col='red')

    #return()
}


#'
#' Fit 2 parameter model with a range of different hurst exponents, and use the 'best'
#' Notice how the 'best' hurst exponent will be model-dependent
fit_multiple_hurst <-function(myrast, hurst_values = c(0.3, 0.4, 0.5+(3.6-0.5)*seq(0,1,len=10)**1.5)){

    # Fit myrast with a range of hurst exponents
    #X = mcmapply(fit_slip_parameters, 
    #    m1 = replicate(length(hurst_values), as.matrix(myrast), simplify=FALSE), 
    #    local_hurst = as.list(hurst_values), 
    #    mc.cores=config_pars$MC_CORES, SIMPLIFY=FALSE, mc.preschedule=FALSE)
    config_pars$check_unique_random_seeds()
    X = clusterMap(cl = config_pars$cluster,
            fun = fit_slip_parameters, 
            m1 = replicate(length(hurst_values), as.matrix(myrast), simplify=FALSE), 
            local_hurst = as.list(hurst_values), 
            .scheduling='dynamic')
    config_pars$check_unique_random_seeds()


    # Extract matrix of fitted pars + values
    fitted_par = matrix(unlist( lapply(X, f<-function(x) x$par)), ncol=2, byrow=T) 
    fitted_values = unlist(lapply(X, f<-function(x) x$value))

    min_ind = which.min(fitted_values)

    if(min_ind%in%c(1,length(hurst_values))){
        # The best fit is an endpoint -- get it
        best_ind = which.min(fitted_values)
        best_hurst = hurst_values[best_ind]
        best_par = fitted_par[best_ind,]
    }else{
        # There should be a smooth relation between the hurst exponent and the
        # goodness-of-fit parameter

        # Solve for the turning point of hurst vs fitted_values
        sfv = splinefun(hurst_values, fitted_values, method='natural')
        out_hurst = try(
            uniroot( f<-function(x) sfv(x, deriv=1), 
                lower=hurst_values[min_ind-1], 
                upper=hurst_values[min_ind+1],
                extendInt='no')$root
            )
        if(class(out_hurst)=='try-error'){
            # Occasionally uniroot can fail if fitted_values is not smooth
            out_hurst = hurst_values[min_ind]
        }

        best_hurst = out_hurst

        sf_p1 = splinefun(hurst_values, fitted_par[,1])
        sf_p2 = splinefun(hurst_values, fitted_par[,2])

        best_par = c(sf_p1(best_hurst), sf_p2(best_hurst))
    }

    output = list(fits = X, best_par = best_par, best_hurst = best_hurst)

    return(output)
}

