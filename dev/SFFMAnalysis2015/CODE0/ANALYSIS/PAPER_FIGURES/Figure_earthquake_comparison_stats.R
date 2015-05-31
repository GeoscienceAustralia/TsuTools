source('load_images.R')


# Function to compute quantile confidence intervals based on binomial theorem
quantile_CI<-function(x, Quant=0.5, level=0.95){
    alpha = (1-level) 
    n = length(x) # Size of data
    allQuants = pbinom(0:n, n, Quant) # All quantiles
    lower_index = max(which(allQuants <= alpha/2))
    upper_index = min(which(allQuants >= (1-alpha/2)))
    if(lower_index == -Inf | upper_index == n+1){
        message=paste('Cannot compute Quantile', 
                    Quant,' at ', level*100, '% confidence',
                    ' with ', n, ' data points')
        stop(message) 
    }
    output=list()
    # Interval
    output$ci = sort(x)[c(lower_index,upper_index)]
    # Exact probability (>level)
    output$p = allQuants[upper_index]-allQuants[lower_index]
    return(output)
}
#quantile_CI(50, Quant = 0.5, level = 0.95)


figdir = 'FIGS_FOR_PAPER'
nEqSims = 50
nFFI = 66

# Extract simulated 0.95 and 0.5 quantile for the +80 region of every FFI
ffi_q95_80 = list()
ffi_median_80 = list()
sffm_q95_80 = list()
sffm_median_80 = list()
for(iter in 1:length(envs)){
    attach(envs[[iter]])
    ffi_q95_80[[iter]] = lapply(rod$myrast_stats, f<-function(x) x$quant_95_80_slip)
    ffi_median_80[[iter]] = lapply(rod$myrast_stats, f<-function(x) x$median_80_slip)
    sffm_q95_80[[iter]] = lapply(rod$myrast_simstats, f<-function(x) x$quant_95_80_slip)
    sffm_median_80[[iter]] = lapply(rod$myrast_simstats, f<-function(x) x$median_80_slip)
    detach(envs[[iter]])
}


# Plot here
pdf(paste0(figdir,'/SummaryStatPlots.pdf'),width=12,height=8)
plotExpressions = titles_envs
par(mfrow=c(3,length(envs)))
powTrans = 1.0
par(mar=c(3,0,0,0))
par(oma=c(0,3,1,0.5))

quickCode<-function(med_ffi_80, med_sffm_80, powTrans, RNG=NULL,logPar='xy', addLab=TRUE){
    # Make a workhorse function to avoid copy-paste of code
    LAB_PTS = RNG
    plot(LAB_PTS**powTrans,LAB_PTS**powTrans, col=0,
         axes=FALSE, asp=1, xaxs='i', yaxs='i', 
         xlab="", ylab="", frame.plot=TRUE, log=logPar)
    abline(0, 1, col='red', lty='longdash')
    title(xlab='FFI', ylab='SFFM', line=-1.5, cex.lab=1.6)
    axis(side=1, at=LAB_PTS**powTrans, cex.axis=1.2,
                    labels=c("", LAB_PTS[2:(length(LAB_PTS)-1)], "")) 
    if(addLab){
        axis(side=2, at=LAB_PTS**powTrans, 
                    labels=LAB_PTS , 
                    cex.axis=1.2, las=1)
    }else{
        axis(side=2, at=LAB_PTS**powTrans, labels=NA)
    }
    XX = c(t(replicate(ncol(med_sffm_80), med_ffi_80)))
    YY = c(t(med_sffm_80))
    points(XX**powTrans, YY**powTrans, pch='.', col='darkgrey')
    loesSmooth = loess(YY~XX,span=0.7)
    predXX = seq(min(XX), max(XX), len=100)
    predYY = predict(loesSmooth, predXX)
    points(predXX**powTrans, predYY**powTrans, t='l', col='green', lwd=3)
}

quickCode_2<-function(med_ffi_80, med_sffm_80, powTrans, RNG=NULL){
    # Make a workhorse function to avoid copy-paste of code
    if(is.null(RNG)){
        RNG = range(med_ffi_80) #range(c(range(med_ffi_80), range(med_sffm_80)))
    }
    LAB_PTS = pretty(RNG,n=7)
    LAB_PTS2 = seq(-1.5,1.5,len=length(LAB_PTS)) #pretty(c(-3, 3),n=7)
    plot(LAB_PTS**powTrans,LAB_PTS2, col=0,
         axes=FALSE,xaxs='i',yaxs='i',xlab="",ylab="",frame.plot=TRUE)
    abline(h=0,col='red',lty='longdash')
    title(xlab='FFI',line=-1.5,cex.lab=1.5)
    title(ylab='[ FFI - SFFM ] / FFI ',line=-1.5,cex.lab=1.5)
    axis(side=1, at=LAB_PTS**powTrans,labels=LAB_PTS,cex.axis=1.2)
    axis(side=2, las=1) # at=LAB_PTS**powTrans,labels=LAB_PTS,cex.axis=1.2,las=1)

    XX = c(t(replicate(ncol(med_sffm_80), med_ffi_80)))
    YY = (XX-c(t(med_sffm_80)))/XX
    points(XX**powTrans, YY,pch='.')
    loesSmooth = loess(YY~XX,span=0.5)
    predXX = seq(min(XX),max(XX),len=100)
    predYY = predict(loesSmooth,predXX)
    points(predXX**powTrans,predYY,t='l',col='green',lwd=5)
}

min_med_80=list()
max_med_80=list()
med_med_80=list()
quant_raw_med_80=list()
ci_raw_med_80=list()
cv_med_80=list()
rcv_med_80=list()
min_95_80=list()
max_95_80=list()
med_95_80=list()
quant_raw_95_80=list()
ci_raw_95_80=list()
cv_95_80=list()
rcv_95_80=list()
min_skew=list()
max_skew=list()
med_skew=list()
quant_raw_skew=list()
ci_raw_skew=list()
cv_skew=list()
rcv_skew=list()

for(i in 1:length(envs)){
    
    # Plot stats for each model
    par(mfg=c(1,i))
    med_ffi_80 = as.numeric(unlist(ffi_median_80[[i]]))
    med_sffm_80 = matrix(unlist(sffm_median_80[[i]]),ncol= nEqSims, byrow=T)
    #quickCode_2(med_ffi_80,med_sffm_80,powTrans=powTrans, RNG=c(0,30))
    #quickCode(med_ffi_80,med_sffm_80,powTrans=powTrans, RNG=c(0.1 , 1.0, 10, 20, 50), addLab =(i==1))
    quickCode(med_ffi_80,med_sffm_80,powTrans=powTrans, RNG=c(0.1 , 1.0, 10, 100, 140.), addLab =(i==1))
    alg = plotExpressions[[i]]
    mtext(bquote(phi[0.5]^{+80}~.(alg)), cex=1.5, adj=0.05, padj=1.3)
    
    # Summary stats
    min_med_80[[i]] = apply(med_sffm_80, 1, min) # Min synthetic median for each FFI
    max_med_80[[i]] = apply(med_sffm_80, 1, max) #
    med_med_80[[i]] = apply(med_sffm_80, 1, median) #
    quant_raw_med_80[[i]] = quantile(med_sffm_80/matrix(med_ffi_80,ncol=nEqSims,nrow=nFFI),p=c(0.25,0.5,0.75)) # Quantiles of ratio of 'median-significant-slip' stats
    ci_raw_med_80[[i]] = quantile_CI(med_sffm_80/matrix(med_ffi_80,ncol=nEqSims,nrow=nFFI)) 
    cv_med_80[[i]] = apply(med_sffm_80, 1, f<-function(x) sd(x)/mean(x))
    rcv_med_80[[i]] = apply(med_sffm_80, 1, f<-function(x) IQR(x)/median(x))

    par(mfg=c(2,i))
    q95_ffi_80 = as.numeric(unlist(ffi_q95_80[[i]]))
    q95_sffm_80 = matrix(unlist(sffm_q95_80[[i]]),ncol=nEqSims,byrow=T)
    quickCode(q95_ffi_80,q95_sffm_80,powTrans=powTrans, RNG=c(0.1, 1, 10., 100., 140), addLab=(i==1))
    mtext(bquote(phi[0.95]^{+80}~.(alg)), cex=1.5, adj=0.05, padj=1.3)
    
    # Summary stats
    min_95_80[[i]] = apply(q95_sffm_80, 1, min) # Min synthetic q95 for each FFI
    max_95_80[[i]] = apply(q95_sffm_80, 1, max) # 
    med_95_80[[i]] = apply(q95_sffm_80, 1, median)
    quant_raw_95_80[[i]] = quantile(q95_sffm_80/matrix(q95_ffi_80,ncol=nEqSims,nrow=nFFI),p=c(0.25,0.5,0.75))# Quantiles of ratio of 'high-significant-slip' stats
    ci_raw_95_80[[i]] = quantile_CI(q95_sffm_80/matrix(q95_ffi_80,ncol=nEqSims,nrow=nFFI))
    cv_95_80[[i]] = apply(q95_sffm_80, 1, f<-function(x) sd(x)/mean(x))
    rcv_95_80[[i]] = apply(q95_sffm_80, 1, f<-function(x) IQR(x)/median(x))

    par(mfg=c(3,i))
    skew_ffi_80 = med_ffi_80/q95_ffi_80
    skew_sffm_80 = med_sffm_80/q95_sffm_80
    quickCode(skew_ffi_80,skew_sffm_80,powTrans=1, RNG=c(0.1,0.2, 0.4, 0.6, 0.8, 1.0),log='', addLab=(i==1))
    mtext(bquote(chi[0.5/0.95]^{+80} ~.(alg)), cex=1.4, adj=0.05, padj=1.3)

    # Summary stats
    min_skew[[i]] = apply(skew_sffm_80, 1, min)
    max_skew[[i]] = apply(skew_sffm_80, 1, max)
    med_skew[[i]] = apply(skew_sffm_80, 1, median)
    quant_raw_skew[[i]] = quantile(skew_sffm_80/matrix(skew_ffi_80, ncol=nEqSims,nrow=nFFI),p=c(0.25,0.5,0.75))# Quantiles of ratio of 'flatness-significant-slip' stats
    ci_raw_skew[[i]] = quantile_CI(skew_sffm_80/matrix(skew_ffi_80, ncol=nEqSims,nrow=nFFI))

    cv_skew[[i]] = apply(skew_sffm_80, 1, f<-function(x) sd(x)/mean(x))
    rcv_skew[[i]] = apply(skew_sffm_80, 1, f<-function(x) IQR(x)/median(x))
}
dev.off()

# Check how often we 'envelope' the stat
envelope_med = rep(NA,length(envs))
envelope_95 = rep(NA,length(envs))
envelope_skew = rep(NA,length(envs))
for(i in 1:length(envs)){
    # Enveloping of median
    med_ffi_80 = as.numeric(unlist(ffi_median_80[[i]]))
    envelope_med[i] = mean((max_med_80[[i]]> med_ffi_80)*(min_med_80[[i]]<med_ffi_80))*100
    
    # 95
    q95_ffi_80 = as.numeric(unlist(ffi_q95_80[[i]]))
    envelope_95[i] = mean((max_95_80[[i]]>q95_ffi_80)*(min_95_80[[i]]<q95_ffi_80))*100.

    # Skew
    skew_ffi_80 = med_ffi_80/q95_ffi_80
    envelope_skew[i] = mean( (max_skew[[i]]>skew_ffi_80)*(min_skew[[i]]<skew_ffi_80) )*100.
}
envelope_med = round(envelope_med)
envelope_95 = round(envelope_95)
envelope_skew = round(envelope_skew)

## Convert to a latex table
myAlg = titles_latex #c('$S_{C}$', '$S_{A}$', '$S_{AF}$', '$S_{CF}$')
myTab = data.frame(Algorithm=myAlg, Envelope.med=envelope_med, Envelope.q95=envelope_95, Envelope.skew=envelope_skew, stringsAsFactors=FALSE)
write.table(myTab,file=paste0(figdir,'/myTab2'),sep=" & ",row.names=FALSE, quote=FALSE, eol='\\\\ \n')

## How well do the models fit [consider ratio of median of SFFMs for each FFI, to FFI value]
for(i in 1:length(envs)) print(summary(med_med_80[[i]]/med_ffi_80))
for(i in 1:length(envs)) print(summary(med_95_80[[i]]/q95_ffi_80))
for(i in 1:length(envs)) print(summary(med_skew[[i]]/skew_ffi_80))

# Convert to stats -- focus on median of the SFFM vs FFI statistics, and the quantiles of the latter
# We compute statistics on ALL INDIVIDUAL SIMULATIONS vs FFIs, despite other explorations above
allQuants = array(NA,dim=c(length(envs), 3, 3))
#dimnames(allQuants) = list(names_envs, c('MedianSigSlip','q95SigSlip','SkewSigSlip'), c('q0.25', 'q0.5', 'q0.75'))
dimnames(allQuants) = list(names_envs, c('MedianSigSlip','q95SigSlip','SkewSigSlip'), c('lower', 'q0.5', 'upper'))

# First dimension is k^{-2}_{*} model, second is the comparison statistic, third holds the quantiles
for(i in 1:length(envs)){
    allQuants[i,1,] = quant_raw_med_80[[i]] #quantile(med_med_80[[i]]/med_ffi_80,p=c(0.25,0.5,0.75))
    allQuants[i,2,] = quant_raw_95_80[[i]]#quantile(med_95_80[[i]]/q95_ffi_80,p=c(0.25,0.5,0.75))
    allQuants[i,3,] = quant_raw_skew[[i]] #quantile(med_skew[[i]]/skew_ffi_80,p=c(0.25,0.5,0.75))

    #allQuants[i,1,] = c(ci_raw_med_80[[i]]$ci[1], quant_raw_med_80[[i]][2], ci_raw_med_80[[i]]$ci[2]) 
    #allQuants[i,2,] = c(ci_raw_95_80[[i]]$ci[1],  quant_raw_95_80[[i]][2],  ci_raw_95_80[[i]]$ci[2]) 
    #allQuants[i,3,] = c(ci_raw_skew[[i]]$ci[1],   quant_raw_skew[[i]][2],   ci_raw_skew[[i]]$ci[2]) 
}

# Format the summary statistics for printing
allQuantsR = format(round(allQuants,2),nsmall=2)
p1 = array(' (',dim=dim(allQuants[,,2]))
p2 = array(',',dim=dim(allQuants[,,2]))
p3 = array(') ',dim=dim(allQuants[,,2]))
allQuantsSum = paste(allQuantsR[,,2],p1, allQuantsR[,,1], p2, allQuantsR[,,3],p3,sep="")
allQuantsSummary = matrix(allQuantsSum, ncol=length(envs), byrow=T)
colnames(allQuantsSummary) = myAlg
rownames(allQuantsSummary) = c('MedianSigSlip','q95SigSlip','SkewSigSlip')
write.table(t(allQuantsSummary),file=paste0(figdir,'/myTab3'), sep=" & ", eol='\\\\ \n',quote=FALSE)

