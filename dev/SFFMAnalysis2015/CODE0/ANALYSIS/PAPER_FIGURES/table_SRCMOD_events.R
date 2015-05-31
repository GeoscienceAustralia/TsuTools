
######################################################################################
#
# Make a table with the information on FFI from SRCMOD that we used
# 
# Had to go back into the data to extract the SRCMOD info flag, which was output to the following table


eventInfoTable = read.table('../../DATA/SRCMOD_PROCESSED/SRCMOD_SUBDUCTION_DATA/SrcModTable_withTags.csv',
    sep=",", header=T)

## Make a column which records how we used the data
# 1 = used for kcx/kcy estimation + slip comparison, 2 = also used for kcf
# tsunami simulation, 3 = one of 8 events used for all k* models

useFlag = rep(1,66) # All events used for earthquake stages

# Find events used for 62 tsunami simulations

tsuSim = basename(Sys.glob('../SLIP_MODELLING/OUTPUTS_2/DEFORMATION/*S_GCF*/S_*'))

tsuSim = sapply(tsuSim, f<-function(x) as.numeric(strsplit(x, '_')[[1]][5]))
useFlag[tsuSim] = 2

# Find events used for 16 tsunami simulations

tsu8Sim = basename(Sys.glob('../SLIP_MODELLING/OUTPUTS_2/DEFORMATION_COPY*/142769980501406_S_GA__NST_abs_SSD_none_RCS_TRUE_hurst_1_noise_gaussian/S_*'))
tsu8Sim = sapply(tsu8Sim, f<-function(x) as.numeric(strsplit(x, '_')[[1]][5]))
useFlag[tsu8Sim] = 3

outTab=paste(eventInfoTable$mytag, eventInfoTable$Location, eventInfoTable$Data, 
             eventInfoTable$Mw, eventInfoTable$scrmodTag, useFlag,
             sep=" & ")

cat(outTab, file=paste0('FIGS_FOR_PAPER/mytab5'), sep=" \\\\ \n")
