source('load_images.R')


figdir='FIGS_FOR_PAPER'

#
plotExpressions = titles_envs # From load_images.R 
allCols = c('blue', 'brown', 'orange', 'purple', 'yellow', 'pink', 'red', 'green')
# Try to set line dots/dashes so they look unique even in greyscale
allLTYs = c('111111', '222222', '333333', '444444', '21212121', '31313131', '41414141', '42424242')
for(i in 1:length(envs)){
    pdf(paste0(figdir,'/kc_vs_Mw_', names_envs[i],'.pdf'),width=10,height=5)
    par(mfrow=c(1,2))
    plot(envs[[i]]$rod$myrasts_Mw, log10(envs[[i]]$rod$p1$kcx), xlab="",ylab="",las=1)
    abline(coef(envs[[i]]$rod$lmx),col='black',lwd=5) 
    for(j in 1:length(envs)){
        if(j==i) next
        abline(coef(envs[[j]]$rod$lmx),col=allCols[j],lty=allLTYs[j],lwd=2) 
    }
    title(xlab=bquote(M[w]),cex.lab=1.3)
    title(ylab=bquote(log[10](k[cx])),cex.lab=1.3,line=2.5)
    grid()
    
    plot(envs[[i]]$rod$myrasts_Mw, log10(envs[[i]]$rod$p1$kcy), xlab="",ylab="",las=1)
    abline(coef(envs[[i]]$rod$lmy),col='black',lwd=5 )
    for(j in 1:length(envs)){
        if(j==i) next
        abline(coef(envs[[j]]$rod$lmy),col=allCols[j],lty=allLTYs[j],lwd=2) 
    }
    title(xlab=bquote(M[w]),cex.lab=1.3)
    title(ylab=bquote(log[10](k[cy])),cex.lab=1.3,line=2.5)
    grid()

    # Add a legend, with the fitted curve for envs[[i]] in black solid line
    legendCols = allCols
    legendCols[i] = 'black'
    legendLWD = rep(2,length(envs))
    legendLWD[i] = 5
    legendLTY = allLTYs
    legendLTY[i] = 'solid'
    legend('bottomleft', 
           legend=plotExpressions, 
           col=legendCols, lwd=legendLWD, lty=legendLTY,bg='white')

    alg=plotExpressions[[i]]
    title(main=bquote("Best-fit corner wavenumber vs Mw for algorithm " ~ .(alg)),cex.main=1.5,outer=TRUE,line=-2)

    dev.off()

}

