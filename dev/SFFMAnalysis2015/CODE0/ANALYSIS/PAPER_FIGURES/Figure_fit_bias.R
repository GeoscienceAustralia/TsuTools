source('load_images.R')
source('../SLIP_MODELLING/calc_kcx_kcy_1D.R', chdir=TRUE)

# Extract the fitted parameters and their 'repeated-simulate-fit' experiments
allFitsStore = list()
allFitsTestStore = list()
allFitsTestSeed = list()

len_par = 2 # kcx, kcy
nsim = 10 # Number of synthetic test fit cases
nEq = 66 # Number of FFI earthquakes

for(iter in 1:length(model_images)){
    attach(envs[[iter]])

    allFitsStore[[iter]] = matrix(unlist(lapply(rod$all_fits, f<-function(x) x$par)), ncol=2, byrow=T)
    allFitsTestStore[[iter]] = array(NA, dim=c(nEq, len_par, nsim))
    allFitsTestSeed[[iter]] = array(NA, dim=c(nEq, nsim))

    for(j in 1:nEq){
        allFitsTestStore[[iter]][j, , ] = matrix(unlist(lapply(testAllFit[[j]], f<-function(x) x$par)), ncol=2)
        allFitsTestSeed[[iter]][j, ] = unlist(lapply(testAllFit[[j]], f<-function(x) paste(attr(x, 'initial_seed'),collapse="") ))
    }

    detach(envs[[iter]])
}

par(mfrow=c(2,2))
resX = unlist(lapply(envs[[1]]$rod$myrasts, f<-function(x) res(x)[1]))
resY = unlist(lapply(envs[[1]]$rod$myrasts, f<-function(x) res(x)[2]))
NX = envs[[1]]$rod$p1$rastNX
NY = envs[[1]]$rod$p1$rastNY

# Get the '1D fit' values for reference
rastsToFit = list()
rastPars = list()
fitted_1D_pars = list()
for(iter in 1:length(model_images)){
    attach(envs[[iter]]) 
    
    # Take the first synthetic raster for each inversion
    rastsToFit[[iter]] = lapply(rod$sim_rasts, f<-function(x) x[[1]])

    # Record its true parameters, used to compute it
    rastPars[[iter]] = matrix(unlist(lapply(rod$sim_rast_Pars, f<-function(x) x[1,])), ncol=2,byrow=T)

    # Fit the raster with the 1D method
    fitted_1D_pars[[iter]] = matrix(unlist( lapply(rastsToFit[[iter]], f<-function(x) calc_1D_kcx_kcy(x))), ncol=2,byrow=T)
    detach(envs[[iter]])
}

par(mfrow=c(2,4))
RelErr_kcx = list()
RelErr_kcy = list()
RelErr_kcx_1D = list()
RelErr_kcy_1D = list()
options(family='serif')
plotExpressions = titles_envs # Taken from load_images.R

for(i in 1:length(allFitsStore)){

    # Check for uniqueness of synthetic ruptures which were fit. 
    num_unique = length(unique(c(allFitsTestSeed[[i]])))
    num_all = length(c(allFitsTestSeed[[i]]))
    if(num_unique != num_all){
        print(paste0('Case ', i, ' is not unique'))
        print(paste0('Unique: ', num_unique))
        print(paste0('Total: ', num_all))
    }else{
        print(paste0('Case ', i, ' is GOOD'))
    } 

    figdir = 'FIGS_FOR_PAPER'
    options(scipen=5) # Supress scientific notation
    pdf(paste0(figdir,'/FitTest_', names_envs[i], '.pdf'),width=12,height=5)
    par(mfrow=c(1,2))
    par(mar=c(4,5,3,1))

    ####################################################
    # Make a plot of KCX
    real_kcx_1D = rastPars[[i]][,1]/resX
    fitted_kcx_1D = abs(fitted_1D_pars[[i]][,1]/resX)
    sims_kcxN = allFitsTestStore[[i]][,1,]
    sims_kcx = sims_kcxN/resX
    real_kcxN = allFitsStore[[i]][,1]
    real_kcx = real_kcxN/resX
    plotRange = c(0.0005, 0.05)

    plot(plotRange, plotRange,log='xy', col=0, xlab="", ylab="", las=1)
    title(ylab=bquote('Fitted' ~ k[cx]), line=3.5, cex.lab=1.3)
    title(xlab=bquote('True'~ k[cx]), line=3, cex.lab=1.3)
    alg = plotExpressions[[i]]
    title(main = bquote("Fitting"~k[cx]~ "with algorithm" ~ .(alg)), line=1, cex.main=1.5)
    points(replicate(nsim, real_kcx), sims_kcx, cex=0.3, pch=19)
    grid(nx=NULL, equilogs=T, lty='dashed')
    abline(0, 1, col=2)


    # Get the 'values we would have fitted' with the 1D approach
    points(real_kcx_1D, fitted_kcx_1D, col='green', pch=3)

    # Error
    RelErr_kcx[[i]] = ((sims_kcx - replicate(nsim, real_kcx))/replicate(nsim,real_kcx))
    RelErr_kcx_1D[[i]] = ((fitted_kcx_1D - real_kcx_1D)/real_kcx_1D)
    legend('topleft', c('Stochastic Optimization', 'Direct Fit to Eqn 3', 'y=x'), pch=c(19, 3, NA), 
           col=c('black', 'green', 'red'), lty=c(NA,NA,1),bg='white',cex=1.3)
   
    ####################################################### 
    # Make a plot of KCY [basically a paste of code above]
    real_kcy_1D = rastPars[[i]][,2]/resY
    fitted_kcy_1D = abs(fitted_1D_pars[[i]][,2]/resY)
    sims_kcyN = allFitsTestStore[[i]][,2,]
    sims_kcy = sims_kcyN/resX
    real_kcyN = allFitsStore[[i]][,2]
    real_kcy = real_kcyN/resX
    plot(plotRange, plotRange, log='xy', col=0, xlab="", ylab="", las=1)
    title(ylab=bquote('Fitted' ~ k[cy]), line=3.5, cex.lab=1.3)
    title(xlab=bquote('True'~ k[cy]), line=3,cex.lab=1.3)
    title(main=bquote("Fitting "~k[cy]~ "with algorithm" ~ .(alg)), line=1, cex.main=1.5)
    points(replicate(nsim, real_kcy), sims_kcy, cex=0.3, pch=19)
    grid(nx=NULL, equilogs=T, lty='dashed')
    abline(0, 1, col=2)

    # Get the 'values we would have fitted' with the 1D approach
    points(real_kcy_1D, fitted_kcy_1D, col='green', pch=3)
    legend('topleft', c('Stochastic Optimization', 'Direct Fit to Eqn 3', 'y=x'), pch=c(19, 3, NA), 
           col=c('black', 'green', 'red'), lty=c(NA,NA,1), bg='white', cex=1.3)

    # Errors
    RelErr_kcy[[i]] = ((sims_kcy - replicate(nsim, real_kcy))/replicate(nsim,real_kcy))
    RelErr_kcy_1D[[i]] = ((fitted_kcy_1D - real_kcy_1D)/real_kcy_1D)
    dev.off()
}

# Here we compute the median relative errors of the 1D and stochastic methods
# Shows how far off from the true value (relatively) we usually are
kcx_rel_err = lapply(RelErr_kcx, f<-function(x) median(abs(x)))
kcy_rel_err = lapply(RelErr_kcy, f<-function(x) median(abs(x)))
kcx_1d_rel_err = lapply(RelErr_kcx_1D, f<-function(x) median(abs(x)))
kcy_1d_rel_err = lapply(RelErr_kcy_1D, f<-function(x) median(abs(x)))

# Here we compute the median relative errors of the 1D and stochastic methods
# Shows bias
kcx_rel_bias = lapply(RelErr_kcx, f<-function(x) median(x))
kcy_rel_bias = lapply(RelErr_kcy, f<-function(x) median(x))
kcx_1d_rel_bias = lapply(RelErr_kcx_1D, f<-function(x) median(x))
kcy_1d_rel_bias = lapply(RelErr_kcy_1D, f<-function(x) median(x))


# Look at the frequency of 'bad' fits where the relative absolute error is >
# 50%
kcx_rel_err_gt50pc = lapply(RelErr_kcx, f<-function(x) 1 - mean(abs(x)<0.5))
kcx_1d_rel_err_gt50pc = lapply(RelErr_kcx_1D, f<-function(x) 1 - mean(abs(x)<0.5))
kcy_rel_err_gt50pc = lapply(RelErr_kcy, f<-function(x) 1 - mean(abs(x)<0.5))
kcy_1d_rel_err_gt50pc = lapply(RelErr_kcy_1D, f<-function(x) 1 - mean(abs(x)<0.5))



fit_stats_output = data.frame(
                    kcx_relative_error = unlist(kcx_rel_err),
                    kcy_relative_error = unlist(kcy_rel_err),
                    kcx_relative_bias = unlist(kcx_rel_bias),
                    kcy_relative_bias = unlist(kcy_rel_bias),
                    kcx_1d_relative_error = unlist(kcx_1d_rel_err),
                    kcy_1d_relative_error = unlist(kcy_1d_rel_err),
                    kcx_1d_relative_bias = unlist(kcx_1d_rel_bias),
                    kcy_1d_relative_bias = unlist(kcy_1d_rel_bias))

rownames(fit_stats_output) = names_envs
write.table(fit_stats_output, file='fitting_method_errors.csv', sep=",", col.names=NA)

#####################################################################
# Compare the fitted corner wavenumbers with the wavenumbers which can
# be extracted from the rasters
x_wavenumber_lower = 1/(NX*resX)
x_wavenumber_upper = floor(NX/2)/(NX*resX)
y_wavenumber_lower = 1/(NY*resY)
y_wavenumber_upper = floor(NY/2)/(NY*resY)


fitted_x_wavenumbers = list()
fitted_y_wavenumbers = list()
for(i in 1:length(envs)){
    fitted_x_wavenumbers[[i]] = allFitsStore[[i]][,1]/resX
    fitted_y_wavenumbers[[i]] = allFitsStore[[i]][,2]/resY
}

# Fraction of wavenumbers below FFI lower limit, and above FFI upper limit
for(i in 1:8) print(mean(fitted_x_wavenumbers[[i]] <= x_wavenumber_lower*0.5))
for(i in 1:8) print(mean(fitted_x_wavenumbers[[i]] <= x_wavenumber_lower))
for(i in 1:8) print(mean(fitted_x_wavenumbers[[i]] > x_wavenumber_lower))
for(i in 1:8) print(mean(fitted_x_wavenumbers[[i]] > x_wavenumber_upper))

for(i in 1:8) print(mean(fitted_y_wavenumbers[[i]] <= y_wavenumber_lower*0.5))
for(i in 1:8) print(mean(fitted_y_wavenumbers[[i]] <= y_wavenumber_lower))
for(i in 1:8) print(mean(fitted_y_wavenumbers[[i]] > y_wavenumber_lower))
for(i in 1:8) print(mean(fitted_y_wavenumbers[[i]] > y_wavenumber_upper))


# Fraction of wavenumbers < 4x lower FFI limit
for(i in 1:8) print(mean(fitted_x_wavenumbers[[i]] <= 4*x_wavenumber_lower))
for(i in 1:8) print(mean(fitted_y_wavenumbers[[i]] <= 4*y_wavenumber_lower))

## 
RelErr_kcx_rng = lapply(RelErr_kcx, f<-function(x) apply(x, 1, f<-function(x) diff(range(x))))
RelErr_kcy_rng = lapply(RelErr_kcy, f<-function(x) apply(x, 1, f<-function(x) diff(range(x))))

RelErr_kcx_max = lapply(RelErr_kcx, f<-function(x) apply(x, 1, f<-function(x) max(x)))
RelErr_kcy_max = lapply(RelErr_kcy, f<-function(x) apply(x, 1, f<-function(x) max(x)))

