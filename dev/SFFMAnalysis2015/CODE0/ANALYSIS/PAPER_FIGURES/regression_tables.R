source('load_images.R')
# Extract data
figdir='FIGS_FOR_PAPER'
pdf(paste0(figdir,'/regressionVis.pdf'),width=10,height=10)

coefsX=list()
sigmasX=list()
coefsY=list()
sigmasY=list()
for(i in 1:length(envs)){
    print(' ')
    print('#####')
    print(summary(envs[[i]]$rod$lmx))
    coefsX[[i]]=summary(envs[[i]]$rod$lmx)$coefficients
    sigmasX[[i]]=summary(envs[[i]]$rod$lmx)$sigma
    par(mfrow=c(2,2))
    plot((envs[[i]]$rod$lmx))
    print(' ')
    print(summary(envs[[i]]$rod$lmy))
    par(mfrow=c(2,2))
    plot((envs[[i]]$rod$lmy))
    coefsY[[i]]=summary(envs[[i]]$rod$lmy)$coefficients
    sigmasY[[i]]=summary(envs[[i]]$rod$lmy)$sigma
    print(' ')
}
dev.off()


# Check residuals for normality
for(i in 1:length(envs)){
    print(' ')
    print('####')
    print(names_envs[i])
    print(shapiro.test(envs[[i]]$rod$lmx$resid))
    print(shapiro.test(envs[[i]]$rod$lmy$resid))
}

# Compute correlation of residuals
resid_corr = rep(NA, length(envs))
resid_corr_test = list()
for(i in 1:length(envs)){
    resid_corr[i] = cor(envs[[i]]$rod$lmx$resid, envs[[i]]$rod$lmy$resid)

    resid_corr_test[[i]] = cor.test(envs[[i]]$rod$lmx$resid, envs[[i]]$rod$lmy$resid)

}

## Convert to a latex table
## We adjust the 'titles_latex' defined in 'load_images.R'
myAlg = c()
for(i in 1:length(titles_latex)){
    myAlg = c(myAlg, c(titles_latex[i], '~'))
}

myEqn=rep('', length(myAlg))
myStdErA=rep('',length(myAlg))
myStdErB=rep('',length(myAlg))
resid_corrs = round(c(rbind(resid_corr, resid_corr)), 2)
for(i in 1:length(titles_latex)){
    myBasicEqn='$\\log_{10}(KMODEL) \\sim ACOEFM_{w}+BCOEF+SIGMA\\epsilon$'
    myBasicEqn=gsub("KMODEL", 'k_{cx}', myBasicEqn)
    myBasicEqn=gsub("ACOEF", format(round(coefsX[[i]][2,1],2),nsmall=2), myBasicEqn)
    myBasicEqn=gsub("BCOEF", format(round(coefsX[[i]][1,1],2),nsmall=2), myBasicEqn)
    myBasicEqn=gsub("SIGMA", format(round(sigmasX[[i]],2), nsmall=2), myBasicEqn)
    myEqn[2*i-1] = myBasicEqn
    myStdErA[2*i-1] = format(round(coefsX[[i]][2,2],3), nsmall=3)
    myStdErB[2*i-1] = format(round(coefsX[[i]][1,2],3), nsmall=3)
}
for(i in 1:length(titles_latex)){
    myBasicEqn='$\\log_{10}(KMODEL) \\sim ACOEFM_{w}+BCOEF+SIGMA\\epsilon$'
    myBasicEqn=gsub("KMODEL", 'k_{cy}', myBasicEqn)
    myBasicEqn=gsub("ACOEF", format(round(coefsY[[i]][2,1],2), nsmall=2), myBasicEqn)
    myBasicEqn=gsub("BCOEF", format(round(coefsY[[i]][1,1],2), nsmall=2), myBasicEqn)
    myBasicEqn=gsub("SIGMA", format(round(sigmasY[[i]],2), nsmall=2), myBasicEqn)
    myEqn[2*i] = myBasicEqn
    myStdErA[2*i] = format(round(coefsY[[i]][2,2],3), nsmall=3)
    myStdErB[2*i] = format(round(coefsY[[i]][1,2],3), nsmall=3)
}
myTab=data.frame(Algorithm=myAlg, Fit=myEqn, Gradient.SE=myStdErA, Intercept.SE=myStdErB, 
    resid.Corr = resid_corrs, stringsAsFactors=FALSE)
write.table(myTab,file=paste0(figdir,'/myTab'),sep=" & ",row.names=FALSE, quote=FALSE, eol='\\\\ \n')

