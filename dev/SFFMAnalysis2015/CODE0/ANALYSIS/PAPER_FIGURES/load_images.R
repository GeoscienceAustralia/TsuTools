# Code to load R images + various packages that most figures will need

library(sp)
library(raster)
library(rgdal)
library(stabledist)
library(fitdistrplus)


######################################################
### PLOT 1 -- Shows a real slip surface + 
###           simulated surfaces for each model.

#' Convert the Rimages filename to a shorthand name
model_images_names<-function(model_images_filename){

    # A stamp occurring in the Rimage file name, and corresponding shorthand name
    image_names_match = matrix( c(
        "S_GA__NST_abs_SSD_none_RCS_TRUE_hurst_1_noise_gaussian.Rdata"     , 'abs_G'   , 
        "S_GAF_NST_abs_SSD_gaussian_RCS_TRUE_hurst_1_noise_gaussian.Rdata" , 'absT_G'  ,
        "S_GC__NST_clip_SSD_none_RCS_TRUE_hurst_1_noise_gaussian.Rdata"    , 'clip_G'  ,
        "S_GCF_NST_clip_SSD_gaussian_RCS_TRUE_hurst_1_noise_gaussian.Rdata", 'clipT_G' ,
        "S_SA__NST_abs_SSD_none_RCS_TRUE_hurst_1_noise_stable.Rdata"       , 'abs_S'   ,      
        "S_SAF_NST_abs_SSD_gaussian_RCS_TRUE_hurst_1_noise_stable.Rdata"   , 'absT_S'  ,
        "S_SC__NST_clip_SSD_none_RCS_TRUE_hurst_1_noise_stable.Rdata"      , 'clip_S'  ,   
        "S_SCF_NST_clip_SSD_gaussian_RCS_TRUE_hurst_1_noise_stable.Rdata"  , 'clipT_S'
        ),
        ncol = 2, byrow=TRUE)

    # Search for the stamps in each filename
    output_name = rep(NA, length(model_images_filename))
    for(i in 1:length(image_names_match[,1])){
        match_index = grep(image_names_match[i,1], model_images_filename)
        stopifnot(length(match_index) == 1)
        output_name[match_index] = image_names_match[i,2]
    }

    # Logical checks on the result
    stopifnot(length(output_name) == length(unique(output_name)))
    stopifnot(sum(is.na(output_name)) == 0)

    return(output_name)
}

# Names of the images
model_images = Sys.glob('../SLIP_MODELLING/OUTPUTS_2/Rimages/run_results_*.Rdata')
stopifnot(length(model_images)==8)
# Get associated shorthand names
names_envs = model_images_names(model_images)

# Read all images
envs=list()
for(i in 1:length(model_images)){
    envs[[i]] = new.env()
    print(paste0('Loading ', model_images[i]))
    load(model_images[i], envir=envs[[i]])
    names(envs)[i] = names_envs[i]
}

# Make plot titles for all images
titles_envs = c(
                expression(S['NA-']), expression(S[NAF]),
                expression(S['NC-']), expression(S[NCF]),
                expression(S['SA-']), expression(S[SAF]),
                expression(S['SC-']), expression(S[SCF])
                ) 

titles_latex = c('$S_{NA-}$', '$S_{NAF}$',
                 '$S_{NC-}$', '$S_{NCF}$',
                 '$S_{SA-}$', '$S_{SAF}$',
                 '$S_{SC-}$', '$S_{SCF}$')

