# Figure with an inverted + synthetic slip surfaces


# Get data
source('load_images.R')

# Raster index
ri = 45
figdir = 'FIGS_FOR_PAPER'
dir.create(figdir, showWarnings=FALSE)
#pdf(paste0(figdir,'/Slip_models_Many.pdf'), width=8, height=6)
pdf(paste0(figdir,'/Slip_models.pdf'), width=8, height=6)

#
# Initially I made plots with a range of different seeds, and then chose one for which
# the stable / gaussian models gave a similar qualitative agreement with the FFI, so 
# the readers can focus on the important issues (range of slip, spread of slip, etc),
# rather than whether one-or-other model works better in that particular random case
# (which anyway will not be reflective of the overall behaviour, which is rather shown
#  in the statistics)

#plot_seed = 1:35
plot_seed = 35

for(rSeed in plot_seed){
    #
    # Get 'real' and 'simulated' FFMs for plot
    FFM = envs[[1]]$rod$myrasts[[ri]] # 'REAL' inversion
    sim_FFMs = list() # Hold simulated FFMs
    bfPars = list() # best-fit parameters
    for(iterator in 1:length(envs)){
        attach(envs[[iterator]])
        set.seed(rSeed) # Make plot reproducible + same random phase
        # Numerical-space fitted parameters
        bfPars[[iterator]]=c(rod$p1$kcxN[ri], rod$p1$kcyN[ri])
        # Simulate
        sim_FFMs[[iterator]]=simulate_sffm(bfPars[[iterator]], FFM)
        detach(envs[[iterator]])
    }

    # Get range for plot
    minSim = min(unlist(lapply(sim_FFMs, f<-function(x) min(as.matrix(x)))))
    maxSim = max(unlist(lapply(sim_FFMs, f<-function(x) max(as.matrix(x)))))
    MINP = min(minSim, min(as.matrix(FFM)))
    MAXP = max(maxSim, max(as.matrix(FFM)))

    ##################################################################
    # FIGURE
    par(family='serif')
    par(mfrow=c(3,3))
    par(mar=c(2.8,2.7,2.5,0.1))
    myCol = heat.colors(20)

    image(FFM, zlim=c(MINP, MAXP), col=myCol)
    contour(FFM,add=T,levels=pretty(c(MINP,MAXP)))

    title('Finite Fault Inversion',cex.main=2)
    title(xlab='Along strike (km)',ylab='Down Dip (km)',line=1.8)

    # Choose sub-figure locations
    figLocs = list(c(1,2), c(1,3), c(2,1), c(2,2), c(2,3), c(3,1), c(3, 2), c(3,3), c(3,4))
    titles = titles_envs # Sourced above
    for(i in 1:length(sim_FFMs)){
        par(mfg=figLocs[[i]])
        image(sim_FFMs[[i]], zlim=c(MINP, MAXP), col=myCol)
        contour(sim_FFMs[[i]], add=T, levels=pretty(c(MINP,MAXP)))
        title(main=titles[i], cex.main=2)
        title(xlab='Along strike (km)', ylab='Down Dip (km)', line=1.8)
    }

    ## Add 'legend' type figure
    #par(mfg=c(2,1))
    #plot(c(0,1),c(0,1),col=0,axes=FALSE,frame.plot=FALSE,xaxs='i',yaxs='i')
    #xp = 0.01
    #yp = c(0.95, 0.95-0.3, 0.95-0.45, 0.95-0.6, 0.15)
    #dyp = -0.2
    #labs = list(paste0('Inversion of ', envs[[1]]$rod$myrastMetaData$Location[ri], 
    #                   '\n Mw: ', envs[[1]]$rod$myrastMetaData$Mw[ri], ', ', 
    #                   envs[[1]]$rod$myrastMetaData$Date[ri], '\n ', 
    #                   envs[[1]]$rod$myrastMetaData$Author[ri]))
    #for(i in 1:length(labs)) text(xp, yp[i], label=labs[[i]], offset=0, adj=c(0,1),cex=1.5)

    ## Colorbar
    #library(colorRampPC)
    #myCol_breaks = breaks4image(myCol, seq(MINP,MAXP,len=20), type='linear')
    #plot_colorVec(myCol, breaks=myCol_breaks, 
    #              xleft=0., xright=1., ybottom=0.3, ytop=0.4,
    #              add=TRUE,add_axis=TRUE)
    #text(0.42, 0.12, 'Slip (m)',adj=c(0,1),cex=1.2)
}
dev.off()

