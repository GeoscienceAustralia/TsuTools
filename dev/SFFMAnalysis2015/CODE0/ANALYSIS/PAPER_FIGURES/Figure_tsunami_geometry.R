
####################################################################################
####################################################################################
#
# Plot of tsunami scenario geometry profile
#
#

rm(list=ls(all=T))
source('tsunami_depth_profile.R')


figdir='FIGS_FOR_PAPER'
options(scipen=5) # Supress scientific notation
pdf(paste0(figdir,'/TsunamiTransect.pdf'),width=7,height=5)
par(mar=c(4,4,1,2))
x=seq(-300000, 300000, by=1000)
y=depthFun(x)


minPloty=-20000
maxPloty=5000*0.01
plot(x, -y, ylim=c(minPloty,maxPloty),xlab='Distance (m)', ylab='Elevation (m MSL)',xaxs='i',yaxs='i')
polygon(c(x,min(x)), c(-y, max(-y)),col=rgb(t(col2rgb('skyblue')), alpha=100, maxColorValue=255))
polygon(c(x,max(x), min(x)), c(-y, minPloty, minPloty),col=rgb(t(col2rgb('brown')), alpha=100, maxColorValue=255))
points(c(0,300000), c(-6000, -6000), t='l',lty='dotted')
points(c(0,300000), c(-6000, -6000 - atan(10/180*pi)*300000), t='l',lty='dashed')
text(15000, -6600, expression(delta),cex=1.4)
#theta=seq(0,10/180*pi,len=200)
#r=10000
#points(r*cos(theta),-6000.-r*sin(theta),t='l')
text(-220000, -2000, 'Ocean',cex=2)
text(-220000, -8000, 'Earth',cex=2)
L=30000
Htop=L*(0.04+atan(10/180*pi))
rupTop=-6000+L*0.04-Htop
arrows(L, -6000+L*0.04, L, rupTop, code=3,angle=90,length=0.05,lwd=2) 
text(L+40000, -6000-3000, expression(h[top]),cex=2)

arrows(L, rupTop, L+50000, rupTop-atan(10/180*pi)*50000, col='red',lwd=2,code=2)
text(L+30000, rupTop-5000, 'Rupture Plane',cex=2,adj=c(0,0))
dev.off()

