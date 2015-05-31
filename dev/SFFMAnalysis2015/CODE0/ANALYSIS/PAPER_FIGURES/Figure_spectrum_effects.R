# Plot of the spectrum which illustrates the issues with clipping
source('load_images.R')

#
# BEWARE: To simulate an unclipped model, I have to make up config_pars
#
config_pars = list()
config_pars$negative_slip_removal_method = 'none'
config_pars$spatial_slip_decay = 'none'
config_pars$RECENTRE_SLIP = FALSE
config_pars$hurst = 1.0
config_pars$noise_distribution = 'gaussian'
config_pars$model_fourierfun = function(kx, ky, reg_par, hurst = config_pars$hurst){
    return( 
        (1.0 + ((kx/reg_par[1])**2 + (ky/reg_par[2])**2)**(1.0+hurst))**(-0.5)
    )
}

source('../SLIP_MODELLING/high_level_analysis.R', chdir=TRUE)
source('../SLIP_MODELLING/calc_kcx_kcy_1D.R', chdir=TRUE)



figdir = 'FIGS_FOR_PAPER'
pdf(paste0(figdir,'/ClippedSpectrum.pdf'), width=8, height=6)
set.seed(1)
nx = 20
ny = 12
template_event = raster(matrix(runif(nx*ny),ncol=nx), xmn=0, xmx=200, ymn=0, ymx=120)

# Force to have mean = 1
mte = mean(as.matrix(template_event))
template_event = template_event/mte

# Set kcx/kcy
dx = res(template_event)[1]
dy = res(template_event)[2]
kcx = 1/100
kcxN = kcx*dx
kcy = 1/100
kcyN = kcy*dy

# Make random matrices -- both should have mean of 1
random_unclipped = simulate_sffm(c(kcxN,kcyN), tg_mat=template_event)
random_clipped = random_unclipped*(random_unclipped>0)
random_clipped = random_clipped/mean(as.matrix(random_clipped))


# Numerical wavenumber info
tmp = get_wavenumbers(random_clipped)
kxN = tmp[[1]]
kyN = tmp[[2]]

# kr = physical space
kr = (1+ ((kxN/(kcxN))**2 + (kyN/(kcyN))**2)**2)**0.5
kUseful = ((kxN/dx)**2+(kyN/dy)**2)**0.5

## "Least-squares fit to clipped model" kr
source('../SLIP_MODELLING/calc_kcx_kcy_1D.R', chdir = TRUE)
kBiasN = calc_1D_kcx_kcy(random_clipped)
krBias = (1+((kxN/kBiasN[1])**2 + (kyN/kBiasN[2])**2)**2)**0.5

kUnbiasN = calc_1D_kcx_kcy(random_unclipped)

# Print out fitted values
print('The following values fit to the unclipped model should be very close to 1/100')
print(kUnbiasN/c(dx,dy))
print('The following values fit to the clipped model should show strong bias (the correct answer is 1/100)')
print(kBiasN/c(dx,dy))

#par(mfrow=c(2,2))
par(family='serif')
layout(matrix(c(1,2,3,3,4,4), ncol=2, byrow=T), 
       width=c(0.5,0.5), height=c(0.5,0.05,0.5))
par(mar=c(4,4,3,1))
RNGE = range(c(as.matrix(random_unclipped), as.matrix(random_clipped)))
myCol = heat.colors(20)
library(colorRampPC)
myColBreaks = breaks4image(myCol, RNGE, type='linear')

image(random_unclipped, zlim=RNGE, 
      main='Unclipped Model Slip' ,
      col=myCol, breaks=myColBreaks,cex.main=1.5,xlab="",ylab="")
contour(random_unclipped, levels=pretty(RNGE, n=6),add=T)
legend('topright', 'A)', bg='white', cex=2, adj=c(0.5,0.5))
title(xlab='Along strike (km)', ylab='Down dip (km)', cex=1.5, line=1.8)

image(random_clipped, zlim=RNGE, 
      main='Clipped Model Slip', 
      col=myCol, breaks=myColBreaks, 
      cex.main=1.5, xlab="", ylab="")
title(xlab='Along strike (km)', ylab='Down dip (km)', 
    cex=1.5, line=1.8)
contour(random_clipped, levels=pretty(RNGE,n=6), add=T)
legend('topright', 'B)', bg='white', cex=2, adj=c(0.5,0.5))

par(mar=c(1,1,0,1))
plot(c(0,1),c(0,1), frame.plot=FALSE, 
    axes=FALSE, col=0, xlab="", ylab="")
plot_colorVec(myCol, breaks=myColBreaks, add=T, add_axis=TRUE, 
    labels=pretty(RNGE,n=6))
mtext('Slip (m)', side=1, line=1, adj=0.49)

ModUnclipped = Mod(fft(as.matrix(random_unclipped)))
ModClipped = Mod(fft(as.matrix(random_clipped)))

par(mar=c(4.5,4,2,1))
plot(kUseful, ModUnclipped, log='y', col='red',
    ylim=c(0.1, max(ModClipped)), xlab="", ylab='',
    t='h', las=1)
points(kUseful, ModClipped,t='p',pch=19)

sort_kr = sort(kr,index.return=T)
sort_ku = kUseful[sort_kr$ix]
EQN = 1*nx*ny/sort_kr$x
points(sort_ku,EQN,t='l',col='blue',lty='dashed')

mtext(expression(sqrt(kx^2 + ky^2)),side=1, line=3)
title(ylab='DFT Amplitude (m)', line=2.6,cex.lab=1.3)
legend('bottomleft', c('Unclipped', 'Clipped', 'Equation 3'), 
       lty=c('solid',NA,'dashed'), col=c('red', 'black', 'blue'), 
       pch=c(NA,19, NA),cex=1.5,bg='white',ncol=2)
legend('topright', 'C)',bg='white',cex=2,adj=c(0.5,0.5))

dev.off()
