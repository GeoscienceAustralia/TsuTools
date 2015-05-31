source('load_images.R')
figdir = 'FIGS_FOR_PAPER'

# Get summary inundation information
inundationSum = read.table('../TSUNAMI_ANALYSIS/AllScenariosSummaryTable.csv',
    sep=",", header=T, stringsAsFactors=FALSE)

# For each synthetic event we want to know it's
# 1) index i1 in envs; 
# 2) index i2 in envs[[i1]]$tsuMod; 
# 3) index i3 in envs[[i1]]$tsuMod[[i2]][[3]]
i1 = rep(NA, length(inundationSum[,1]))
i2 = rep(NA, length(inundationSum[,1]))
i3 = rep(NA, length(inundationSum[,1]))

for(i in 1:length(model_images)){
    ith_image_grep_flag = strsplit(basename(model_images[i]), '_')[[1]][3]
    i1[grep(ith_image_grep_flag, inundationSum[,1])] = i
}

# i2 index comes from matching "the number _i appended to the event name", with
# "the same _i number name thing in the tsuMod object"
#
# Actually we can do it this way
._i_in_eventName = unlist(lapply(basename(as.character(inundationSum[,1])),
                        f<-function(x) as.numeric(strsplit(x,'_')[[1]][5])))

# Get the corresponding integer, as ordered in tsuMod (and rod)
._i_in_tsuModInd = unlist( lapply(basename(envs[[1]]$rod$myrastlist), 
                         f<-function(x) as.numeric(gsub('.tif','', strsplit(x, '_')[[1]][5])) ))

i2 = match(._i_in_eventName, ._i_in_tsuModInd)

# i3 index comes from the OceanInitial_i term
._OceanInitialWord = sapply(as.character(inundationSum[,2]),
                        f<-function(x) strsplit(x,'/',fixed=T)[[1]][5] )
i3 = as.numeric(gsub('OceanInitial_', '', ._OceanInitialWord)) # NA's are FFI

# Extract some stats about the slip. Suggest we get median 80+ and 95th percentile 80+
median_95_80p_slip<-function(rast, signif_thres=0.8, quantileVal=0.95){
    # Utility function to return median+q95 80+ slip, and raw mean slip
    rastM = as.matrix(rast)
    sortedSlip = sort(rastM)
    sortedSlip_CS = cumsum(sortedSlip)
    SignificantSlipVals = sortedSlip[sortedSlip_CS>(1-signif_thres)*max(sortedSlip_CS)]

    return(c(median(SignificantSlipVals), quantile(SignificantSlipVals, p=quantileVal), mean(rastM)))
}

if(FALSE){
    # Quick test that the above function works
    median_95_80p_slip(envs[[3]]$rod$myrasts[[1]])
    # Should have both of these vaues 
    envs[[3]]$rod$myrast_stats[[1]]$median_80_slip
    envs[[3]]$rod$myrast_stats[[1]]$quant_95_80_slip
}

# Now get the statistics for the SFFM and FFI
inundationSlipSum = matrix(NA,nrow=length(inundationSum[,1]), ncol=3)
for(i in 1:length(inundationSum[,1])){
    print(i)
    if(!is.na(i3[i])){
        inundationSlipSum[i,] = median_95_80p_slip( envs[[ i1[i] ]]$tsuMod[[ i2[i] ]][[3]][[i3[i]]] )
    }else{
        # In this case, get it from the FFI raster
        inundationSlipSum[i,] = median_95_80p_slip( envs[[1]]$rod$myrasts[[i2[i]]])
    }
}

## Better -- relate inundation to Mw, then consider the effect of slip
Mw = sapply(basename(as.character(inundationSum[,1])), f<-function(x) as.numeric(gsub('M','', strsplit(x, '_')[[1]][3])))

## Figure out the predicted inundation based on mean slip alone
eventFlag = i2 #i1*100+i2
uniqueEventFlag = unique(eventFlag)
corStore_A = rep(NA,length(uniqueEventFlag))
corStore_B = rep(NA,length(uniqueEventFlag))
corStore_C = rep(NA,length(uniqueEventFlag))
Mw_store = rep(NA,length(uniqueEventFlag))

corMeth = 'spearman'
for(i in 1:length( uniqueEventFlag)){
    keepers = which( (eventFlag==uniqueEventFlag[i]) & (i1==4) ) # Only compute for TSUNAMI_CLIP_TAPER events

    Mw_store[i] = Mw[keepers[1]]
    corStore_A[i] = cor((inundationSum$maxInundation[keepers]), (inundationSlipSum[keepers,1]),method=corMeth,use='pairwise.complete.obs')
    corStore_B[i] = cor((inundationSum$maxInundation[keepers]), (inundationSlipSum[keepers,2]),method=corMeth,use='pairwise.complete.obs')

    corStore_C[i] = cor((inundationSum$maxInundation[keepers]), (inundationSlipSum[keepers,1]/inundationSlipSum[keepers,2]),method=corMeth,use='pairwise.complete.obs')
}
summary(corStore_A)
summary(corStore_B)
summary(corStore_C)

pdf(paste0(figdir,'/MeanSlipvsRunupHeight.pdf'),width=8,height=6)
plot(inundationSlipSum[,3], inundationSum$maxInundation/100.,log='xy',
     col=(is.na(i3)+1)*(i1==4),pch=19, cex=inundationSlipSum[,1]/inundationSlipSum[,2]*2,
     xlab='Mean Slip (m)', ylab='Maximum Runup Height (m)')
legend('bottomright', c('SFFM', 'FFI'), col=c('black', 'red'), pch=19)
abline(0,8,untf=T,col='blue',lty='dashed')
abline(0,1,untf=T,col='blue',lty='dashed')
grid(lty='dashed')
rect(0.008, 25, 0.7, 100,col=rgb(1,1,1,alpha=0.5))
legend(0.08, 50.,  legend=c(0.2, 0.4, 0.6, 0.8), pt.cex=c(0.2,0.4,0.6, 0.8)*2,  title='Flatness Statistic', pch=19,bty='n',ncol=4)
dev.off()


pdf(paste0(figdir,'/MwVsRunupHeight.pdf'),width=10,height=8)
par(mfrow=c(2,2))
par(mar=c(4,4,2,1))
par(family='serif')
plot(Mw, inundationSum$maxInundation/100.,log='xy',
     col=(is.na(i3)+1)*(i1==4),pch=is.na(i3)*19 + is.finite(i3), xlab="",ylab="",
     main=bquote(M[w] ~ 'vs Runup Height'), cex.main=1.8)
title(xlab=bquote(M[w]), ylab='Maximum Runup Height (m)',cex.lab=1.5,line=2.2)

grid(lty='dashed')
legend('bottomright', c('SFFM', 'FFI'), col=c('black', 'red'), pch=c(1,19), bg='white')
#abline(0,8,untf=T,col='blue',lty='dashed')
#abline(0,1,untf=T,col='blue',lty='dashed')

plot(Mw_store, corStore_A,xlab="",ylab="",
     main=bquote( phi[0.5]^{+80} ~ ' and Runup Height'),cex.main=1.8)
title(xlab=bquote(M[w]~'of each SFFM group'),ylab='Spearman Rank Correlation',cex.lab=1.3, line=2.2)
grid()
abline(h=0,col='black')
plot(Mw_store, corStore_B, xlab="",ylab="",
     main=bquote(phi[0.95]^{+80} ~ ' and Runup Height'), cex.main=1.8)
grid()
abline(h=0,col='black')
title(xlab=bquote(M[w]~'of each SFFM group'),ylab='Spearman Rank Correlation',cex.lab=1.3, line=2.2)

plot(Mw_store, corStore_C, xlab="", ylab="",
     main=bquote(chi[0.5/0.95]^{+80} ~ ' and Runup Height'), cex.main=1.8)
title(xlab=bquote(M[w]~'of each SFFM group'),ylab='Spearman Rank Correlation',cex.lab=1.3, line=2.2)
grid()
abline(h=0,col='black')
dev.off()

