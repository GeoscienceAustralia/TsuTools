
# Setup configuration data
configuration_data = list(

   stable_notaper = list(
                        fourierFun = '2parSom', 
                        noise_distribution = 'gaussian',
                        negative_slip_removal_method = 'abs',
                        spatial_slip_decay = 'none',
                        RECENTRE_SLIP=TRUE,
                        hurst=1.0,
                        use_test_data=FALSE)
   
)

# Ensure appropriate parallel random number generator is used (or the numbers
# might not be independent)
RNGkind("L'Ecuyer-CMRG")

# Ensure that names in the configuration_data are unique
configuration_names = names(configuration_data)
stopifnot(length(unique(configuration_names)) == length(configuration_data))

#
config_name = configuration_names[1]
config_pars = list2env(configuration_data[[config_name]])
print(as.list(config_pars))
source('static_configuration_parameters.R', local = config_pars)

source('high_level_analysis.R')
myrasts = lapply(config_pars$myrastlist, f<-function(x) raster(x))

if(FALSE){
    # Fit all the rasters with a range of hurst exponents
    fits = list()
    for(i in 1:length(myrasts)){
        print(i)
        fits[[i]] = fit_multiple_hurst(myrasts[[i]])
    }
}

##################################################
# Test a synthetic case
# This shows how difficult it is to accurately identify the model with varying
# hurst + corner-wavenumbers
test_par = c(0.075, 0.078) # Corner wavenumbers
test_hurst = 1.0

ntests = 50

# Simulate it and fit it
test_fits=list()
test_rasts = list()
for(i in 1:ntests){
    print(i)
    test_rasts[[i]] = simulate_sffm(test_par, myrasts[[1]], 
        local_hurst=test_hurst)
    test_fits[[i]] = fit_multiple_hurst(test_rasts[[i]])
}

fit_hurst = unlist(lapply(test_fits, f<-function(x) x$best_hurst))
# Check if median might be true value
binom.test(sum(fit_hurst<test_hurst), length(fit_hurst), p=0.5)

fit_kcx = unlist(lapply(test_fits, f<-function(x) x$best_par[1]))
# Check if median might be true value
binom.test(sum(fit_kcx<test_par[1]), length(fit_hurst), p=0.5)

fit_kcy = unlist(lapply(test_fits, f<-function(x) x$best_par[2]))
# Check if median might be true value
binom.test(sum(fit_kcy<test_par[2]), length(fit_hurst), p=0.5)

# Strong correlation between fit_kcx and fit_hurst, and strong variation
# This seems to be why it's hard to fit reliably
par(mfrow=c(1,2))
plot(fit_hurst, fit_kcx)
plot(fit_kcx, fit_kcy)


### Some analysis
#
#hursts = unlist(lapply(fits, f<-function(x) x$best_hurst))
#best_par = matrix(unlist(lapply(fits, f<-function(x) x$best_par)), ncol=2)
#delX = unlist(lapply(myrasts, f<-function(x) res(x)[1]))
#delY = unlist(lapply(myrasts, f<-function(x) res(x)[2]))
#
#best_fit = best_par/cbind(delX, delY)
#mw = config_pars$myrastMetaData$Mw
#
#kcx = best_fit[,1]
#kcy = best_fit[,2]
#
#model1 = lm(log10(kcx) ~ mw)
#model2 = lm(log10(kcy) ~ mw)
