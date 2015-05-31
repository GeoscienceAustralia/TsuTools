# Provide functions for basic laplace distribution
library(HyperbolicDist)

dlap<-function(x, scale){
    return(dskewlap(x, Theta=c(scale, scale, 0.)))
}

qlap<-function(p, scale){
    return(qskewlap(p, Theta=c(scale, scale, 0.)))
}

plap<-function(q, scale){
    return(pskewlap(q, Theta=c(scale, scale, 0.)))
}

rlap<-function(n, scale){
    return(rskewlap(n, Theta=c(scale, scale, 0.)))
}
