###############################################################################
#
# configuration parameters which don't change between analyses. These are
# sourced into an environment named 'config_pars' which also contains the
# 'adjustable' configuration parameters
#
# AUTHORS
# Gareth Davies, Geoscience Australia, 2013-2015,
# gareth.davies.ga.code@gmail.com
#
###############################################################################
library(parallel)

# Ensure appropriate parallel random number generator is used (or the numbers
# might not be independent)
RNGkind("L'Ecuyer-CMRG")


###############################################################################
#
# DATA
#
###############################################################################

if(use_test_data == TRUE){
    # This should 'glob' to all the input slip rasters
    rasterList_glob = '../../DATA/TEST_DATA/slipRasts/*.tif' 
    # Metadata for the slip rasters
    rasterMetaDataFile = '../../DATA/TEST_DATA/SrcMod_MetaData.csv'

}else if(use_test_data == FALSE){
    # This should 'glob' to all the input slip rasters
    rasterList_glob = 
        '../../DATA/SRCMOD_PROCESSED/SRCMOD_SUBDUCTION_DATA/slipRasts/*.tif' 
    # Metadata for the slip rasters
    rasterMetaDataFile = 
        '../../DATA/SRCMOD_PROCESSED/SRCMOD_SUBDUCTION_DATA/SrcMod_MetaData.csv'

}else{

    stop('Invalid value for use_test_data')

}


# Define list of rasters
myrastlist = Sys.glob(rasterList_glob)

# Define raster metadata data.frame -- and order the rows to match myrastlist
myrastMetaData = read.csv(rasterMetaDataFile, header=T) # Temporary, is reordered
# The last number in the filename is the row-index of the metadata
rastIndex = unlist(
    lapply(strsplit(basename(myrastlist),'_'), myfun<-function(x) x[5])
)
rastIndex = as.numeric(gsub('.tif', '', rastIndex))
# Final myrastMetaData 
myrastMetaData = myrastMetaData[rastIndex,]

###############################################################################
#
# Prediction of kcx/kcy values from Mw
#
###############################################################################

# Flag for whether we should use the log10(kcx/y)~Mw regression relations to
# simulate kcx/kcy values for synthetic slip models. Otherwise, we just use
# kcx/kcy fitted from the corresponding earthquake.
random_kcx_kcy_simulations = TRUE

# Should we fit log10(kcx) = -0.5Mw + Constant [the self similar model].
# Otherwise, allow the coefficient -0.5 to be estimated from the data.
self_similar_kcx_mw_regression = FALSE

###############################################################################
#
# Computational parameters
#
###############################################################################

# Number of cores to use in parallel
MC_CORES = detectCores()

# Initial set_seed integer -- nice to control random numbers
initial_set_seed_integer = 12345

# Type of cluster ('MPI' on NCI, 'FORK' on my home machine, ...)
cluster_type = 'FORK'

#' Functon to check that each worker on the cluster has a unique random seed
check_unique_random_seeds<-function(verbose=TRUE, print_result=FALSE){

    seed_info = clusterEvalQ(
        cl = cluster, 
        paste0(as.character(.Random.seed), collapse="")
        )

    if(print_result){
        print(seed_info)
    }

    if(length(seed_info) != length(unique(unlist(seed_info)))){
        print(seed_info)
        stop('Error: Non-unique random seeds')    
    }else{
        if(verbose) print('------ CLUSTER RANDOM SEEDS ARE UNIQUE ----------')
    }

}

# Make the cluster if it doesn't exist
if(!exists('cluster')){
    cluster = makeCluster(MC_CORES, cluster_type)

    # Ensure appropriate parallel random number generator is used (or the numbers
    # might not be independent)
    RNGkind("L'Ecuyer-CMRG")

    set.seed(initial_set_seed_integer)

    clusterSetRNGStream(cl = cluster, iseed = .Random.seed)

    # Load packages
    clusterEvalQ(cl = cluster, {
        library(rgdal)
        library(raster)
        library(stabledist)
        library(colorRampPC)
        library(EqSim)
        library(fitdistrplus)
        library(HyperbolicDist)
        library(ismev)
    }
    )

    # Call a random number on each worker
    #.TMP. = clusterMap(cl=cluster, fun=function(x) return(runif(1)), x=1:MC_CORES)
    check_unique_random_seeds()
}


###############################################################################
#
# Other parameters
#
###############################################################################

# lambda for inverse box-cox transformation
box_cox_lambda = 0.2

# Number of synthetic simulations used to investigate earthquake statistics
num_sims_earthquake_statistics = 50 # 50 # must be >= 3

# Number of synthetic tsunami events per FFI
num_synth_tsunami = 10 # 10

# Number of simulations used to test our model fitting algorithm
num_test_fit_sims = 10 # 30

###############################################################################
#
# Output folders
#
###############################################################################
deformation_folder = 'OUTPUTS_2/DEFORMATION'
quickfig_folder = 'OUTPUTS_2/FIG'
pdf_folder = 'OUTPUTS_2/pdfs'
rimages_folder = 'OUTPUTS_2/Rimages'


###############################################################################
#
# VARIOUS SPECTRAL MODELS
#
###############################################################################

switch(fourierFun,
    '2parSom' = {
        # 2 parameter model -- Sommerville
        model_fourierfun<-function(kx, ky, reg_par, hurst = config_pars$hurst){
            return( 
                (1.0 + ((kx/reg_par[1])**2 + (ky/reg_par[2])**2)**(1.0+hurst))**(-0.5)
            )
        }
    },

    #'3parSom' = {
    #    # 3 parameter model -- Sommerville + variable hurst exponent
    #    model_fourierfun<-function(kx, ky, reg_par){
    #        return( 
    #            (1.0 + ((kx/reg_par[1])**2 + (ky/reg_par[2])**2)**(1.0+reg_par[3]))**(-0.5)
    #        )
    #    }
    #},

    stop(paste0(' Function for fourierFun ', fourierFun, ' not defined'))
)

#
# 1 parameter model -- Gallovic
#model_fourierfun<-function(kx,ky,reg_par) (1.0 + ((kx/reg_par[1])**2 + (ky/reg_par[1])**2)**(2))**(-0.5) 

# Another 1 parameter model -- small constant term
#model_fourierfun<-function(kx,ky,reg_par) (1.0e-04 + ((kx/reg_par[1])**2 + (ky/reg_par[2])**2)**(2))**(-0.5) 

# Another 2 parameter model -- variable constant term
#model_fourierfun<-function(kx,ky,reg_par) (reg_par[2] + ((kx/reg_par[1])**2 + (ky/reg_par[1])**2)**(2))**(-0.5) 

# 3 parameter model 
#model_fourierfun<-function(kx,ky,reg_par){
#    tmp=(1.0 + ((kx/reg_par[1])**2 + (ky/reg_par[1])**2)**(2))**(-0.5)
#    tmp=tmp*matrix(rlnorm(length(kx),meanlog=reg_par[2], sdlog=reg_par[3]), ncol=ncol(kx),nrow=nrow(kx))
#    return(tmp)
#}

## Another 3 parameter model
#model_fourierfun<-function(kx,ky,reg_par){
#        kcx=reg_par[1]
#        kcy=reg_par[2]
#        n=reg_par[3]
#        (1.0 + ((kx/kcx)**2 + (ky/kcy)**2)**(-1/n))**(n)
#}

# Yet another 3 parameter model
#model_fourierfun<-function(kx,ky,reg_par) (1.0 + ((kx/reg_par[1])**2 + (ky/reg_par[2])**2)**(-1/reg_par[3]))**(-0.5)
#######################################################

