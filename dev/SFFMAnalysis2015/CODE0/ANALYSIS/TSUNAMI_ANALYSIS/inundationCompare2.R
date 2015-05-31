# R code to compare max inundation lines

read_all_wdl<-function(sourcedir){
    # Read in all the wet-dry lines inside a given sourcedir
    awd=Sys.glob(paste0(sourcedir,'/*/*/*/wet_dry_line.txt'))

    awdL=list()
    for(i in 1:length(awd)){
        awdL[[i]]=read.table(awd[i],colClasses='numeric')
        names(awdL)[i]=awd[i] #strsplit(awd[i],'/')[[1]][4]
    }
    return(awdL)
}

########################################################

quick_plot<-function(sourcedir, figdir){
    # Plot the max inundation lines

    dir.create(figdir,showWarnings=FALSE,recursive=T)

    # Get all wet_dry_line.txt files
    awdL=read_all_wdl(sourcedir)

    # Make plot -- force axis range
    i_max=-Inf
    i_min=Inf
    for(i in 1:length(awdL)){
        i_max=max(i_max, max(awdL[[i]][,1]))
        i_min=min(i_min, min(awdL[[i]][,1]))
    }
    if(i_max-i_min < 20){
        i_min=i_min-10
        i_max=i_max+10
    }
        

    legendName=basename(basename(basename(dirname(dirname(dirname(names(awdL)))))))
    legendCol=1:length(legendName)
    legendLWD=rep(1,length(legendName))
    k=grep('Orig', legendName)
    if(length(k)>0){
        legendCol[k]=1
        legendLWD[k]=3
    }

    png(paste0(figdir,basename(sourcedir),'.png'),width=14,height=7,units='in',res=300)
        layout(matrix(c(1,2),ncol=2), widths=c(5,2))
        plot(c(i_min,i_max), c(0,1000),col='white')
        for(i in 1:length(awdL)){
            if(length(grep('Orig', legendName[i]))==0){
                points(awdL[[i]],t='l',col=legendCol[i],lwd=legendLWD[i])
            }else{
                points(awdL[[i]],t='l',lwd=legendLWD[i],col=legendCol[i])
            }
        }
        plot(c(0,1),col=0,axes=FALSE,xlab="",ylab="")
        legend('center', legendName,
               col=legendCol, lty='solid', lwd=legendLWD, cex=1)
    dev.off()

}

################################################################

summarise_allWDL<-function(allWDL){
    #@ Extract summary stats from the WDL object

    #sources=rep(NA,length(allWDL)*11)
    #events=rep(NA,length(allWDL)*11)
    #minInun=rep(NA,length(allWDL)*11)
    #meanInun=rep(NA,length(allWDL)*11)
    #maxInun=rep(NA,length(allWDL)*11)
    sources=c()
    events=c()
    minInun=c()
    meanInun=c()
    maxInun=c()
    isEnveloped=c()
    upperBuf=c()
    lowerBuf=c()

    for(i in 1:length(allWDL)){
        #sources[(i-1)*11 + 1:11] = names(allWDL)[i]
        #events[(i-1)*11 + 1:11] = names(allWDL[[i]])
        #minInun[(i-1)*11+1:11] = unlist(lapply(allWDL[[i]], f<-function(x) min(x[,1])))
        #meanInun[(i-1)*11+1:11] = unlist(lapply(allWDL[[i]], f<-function(x) mean(x[,1])))
        #maxInun[(i-1)*11+1:11] = unlist(lapply(allWDL[[i]], f<-function(x) max(x[,1])))
        minInun= c(minInun, unlist(lapply(allWDL[[i]], f<-function(x) min(x[,1])), use.names=F))
        meanInun = c(meanInun, unlist(lapply(allWDL[[i]], f<-function(x) mean(x[,1])), use.names=F))
        maxInun = c(maxInun, unlist(lapply(allWDL[[i]], f<-function(x) max(x[,1])), use.names=F))
        sources = c(sources, rep(names(allWDL)[i], length(minInun)-length(sources)))
        events = c(events, names(allWDL[[i]]))

        tmp=is_enveloped(allWDL[[i]])
        isEnveloped=c(isEnveloped, rep(tmp[1], length(minInun)-length(isEnveloped)))
        upperBuf=c(upperBuf,rep(tmp[2],  length(minInun)-length(upperBuf)))
        lowerBuf=c(lowerBuf,rep(tmp[3],  length(minInun)-length(lowerBuf)))

    }

    meanDF=data.frame(Sources=sources, events=events,minInun=minInun,maxInun=maxInun,meanInun=meanInun,
                      isEnveloped=isEnveloped, upperBuf=upperBuf, lowerBuf=lowerBuf, stringsAsFactors=FALSE)
    return(meanDF)
}

#####################################################################

produce_figures_and_data<-function(){

    # Run in these directories
    sourcedirs=Sys.glob('WDOUT/*/*/S_*')

    allWDL=list()
    counter=0
    for(sourcedir in sourcedirs){
        print(sourcedir)
        quick_plot(sourcedir,figdir=paste0('FIG/',basename(dirname(sourcedir)), '/')) 
        counter=counter+1
        allWDL[[counter]]=read_all_wdl(sourcedir)
        names(allWDL)[counter]=sourcedir
    }
    saveRDS(allWDL,file='all_wet_dry_lines.RDS')

}

##################################################################

is_enveloped<-function(event_WDL, returnSynthEnvelope=FALSE, returnMinMax=FALSE){
    # Take a list of wet-dry-lines for an event,
    # and determine whether the 'Orig' event is
    # enveloped by the other ones

    OrigEventInd=grep('OceanInitial_Orig', names(event_WDL))

    # To get determine the enveloping, use approx to get the WD extents at regular intervals
    # 1, 2, 3, ... 999
    re_gridded_eventWDL=list()
    Xout=seq(1,999,by=0.5)
    for(i in 1:length(event_WDL)){
        re_gridded_eventWDL[[i]] = approx(event_WDL[[i]][,2], event_WDL[[i]][,1], 
            xout=Xout, rule=2)$y
        names(re_gridded_eventWDL)[i] = names(event_WDL)[i]

    }

    # Get the 'min' and 'max' output value == max and min inundation
    # respectively (since lower x coordinate = more onshore)
    l = length(re_gridded_eventWDL[[i]])
    minOut = rep(Inf,l)
    maxOut = rep(-Inf,l)
    for(i in 1:length(event_WDL)){
        if(i == OrigEventInd){
            # Skip it
        }else{
            minOut = pmin(minOut,re_gridded_eventWDL[[i]])
            maxOut = pmax(maxOut,re_gridded_eventWDL[[i]])
        }
    }

    if(returnSynthEnvelope){
        return(list(x=Xout,uLim=minOut,lLim=maxOut, orig=re_gridded_eventWDL[[OrigEventInd]]))
    }else if(returnMinMax){
        Upper_max = max(500.-minOut)
        Lower_min = min(500.-maxOut)

        Upper_FFI = max(500-re_gridded_eventWDL[[OrigEventInd]])
        Lower_FFI = min(500-re_gridded_eventWDL[[OrigEventInd]])

        return(c(Upper_max, Lower_min, Upper_FFI, Lower_FFI))

    }else{
        # Get stats on the closest upper and lower boundary-- remember minout is
        # the highest inundation
        upperSpaceBuffer = min((500-minOut)-(500-re_gridded_eventWDL[[OrigEventInd]]))
        lowerSpaceBuffer = -max((500-maxOut)-(500-re_gridded_eventWDL[[OrigEventInd]]))
        if(all(maxOut >= re_gridded_eventWDL[[OrigEventInd]]) & 
           all(minOut <= re_gridded_eventWDL[[OrigEventInd]])){
            return(c(TRUE, upperSpaceBuffer, lowerSpaceBuffer))
        }else{
            return(c(FALSE, upperSpaceBuffer, lowerSpaceBuffer))
        }
    }

}

#########################################################

explore_wd_enveloping_stats<-function(){
    #
    # Check how often the wet-dry lines envelope each other
    #

    # Get summary info, and make sure output is reasonably sorted
    allWDL = readRDS('all_wet_dry_lines.RDS')
    sumAWDL = summarise_allWDL(allWDL)
    sumAWDL = sumAWDL[sort(sumAWDL$Sources,index.return=T)$ix,]
    Orig = grep('_Orig', sumAWDL$events)

    ##################################################################
    # Some plots
    boxplot((500-sumAWDL$meanInun[-Orig])~ sumAWDL$Sources[-Orig],log='y',col='red')
    boxplot((500-sumAWDL$meanInun[Orig]) ~ sumAWDL$Sources[Orig], add=T,col='blue')

    # Compute min/max, and ensure ordering is ok
    #
    rangeGroup = aggregate(sumAWDL$meanInun[-Orig], 
                         by=list(sumAWDL$Sources[-Orig]),
                         f<-function(x) c(min(x), max(x)) )

    if(!all(rangeGroup[[1]]==sumAWDL$Sources[Orig])){
        stop('SORTING ERROR')
    }else{
        print('Sorting seems ok')
    }

    # Extract the max of the mean inundations for each group
    max_sim_mean = 500-rangeGroup[[2]][,1]
    min_sim_mean = 500-rangeGroup[[2]][,2]
    Orig_mean = 500-sumAWDL$meanInun[Orig]
    Orig_max = 500-sumAWDL$minInun[Orig]

    plot(max_sim_mean,log='y',t='o')
    points(min_sim_mean,col=2,t='o')
    points(Orig_mean,t='o',col=3)
    #
    # Note -- rangeGroup[[2]][,1] gives the Smallest x value, corresponding to the MAX inundation
    #
    over_mean = (Orig_mean>max_sim_mean)
    under_mean = (Orig_mean<min_sim_mean)
    # Get earthquake magnitudes
    SourceNames = sumAWDL$Sources[Orig]
    getMw<-function(SourceName) as.numeric(gsub('M', "", strsplit(basename(SourceName), '_')[[1]][3]))
    Mw = sapply(SourceNames, getMw) 

    #####################################
    # EXTRACT COVERAGE STATISTICS 
    # OUTPUT FOR SIXTEEN EVENTS, ALL SFFM
    OrigSum = sumAWDL[Orig,]

    # Names to match different SFFM types
    model_type_matches = c('S_GA_', 'S_GC_', 'S_GAF', 'S_GCF', 'S_SA_', 'S_SC_', 'S_SAF', 'S_SCF')

    EightRunsCover=list()
    for(model_type in model_type_matches){
        inds = grep(model_type, OrigSum[,1], fixed=TRUE)
        inds = intersect(inds, grep('EIGHT_FFI', OrigSum[,1], fixed=TRUE) )

        EightRunsCover[[model_type]] = data.frame(Source=basename(OrigSum[inds,1]),
            isEnveloped=OrigSum[inds,6],
            UpBuf=OrigSum[inds, 7], 
            loBuf=OrigSum[inds, 8])
    }

    ##############################################
    FFI_no = unlist(lapply(strsplit(as.character(EightRunsCover[[1]]$Source), '_'), f<-function(x) x[5]))
    Mw_these = sapply(as.character(EightRunsCover[[1]]$Source),f<-function(x) gsub('M',"", strsplit(x,'_')[[1]][3]) )
    colNames = c('FFI ID', model_type_matches)

    # Export distance between event and envelope
    # For positive errors, we export a positive number, negative for a negative
    # number, zero for enveloped events
    # We also add in the max inundation from the FFI
    outTable = data.frame(FFI = FFI_no, Mw = as.character(Mw_these))
    FFI_inundation = grep(model_type_matches[1], OrigSum[,1], fixed=TRUE)
    FFI_inundation = intersect(FFI_inundation, grep('EIGHT_FFI', OrigSum[,1], fixed=TRUE))
    outTable[['FFI Inundation']] = as.character(round(Orig_max[FFI_inundation]*10,0))
    for(model_type in model_type_matches){
        X = EightRunsCover[[model_type]]
        errs = (-X[,3]*(X[,3]<0) + X[,4]*(X[,4]<0))*10
        reformatInds = (errs != 0)
        errs[reformatInds] = paste('\\textbf{', errs[reformatInds], '}',sep="")
        outTable[[model_type]] = errs
    }

    # Sort with increasing Mw
    newOrder=sort(Mw_these,index.return=T)
    outTable=outTable[newOrder$ix,]

    write.table(outTable,file='EightRunsStats.txt',sep=" & ",
                row.names=FALSE,col.names=TRUE,quote=F,
                eol=' \\\\\n')

    # Some stats which we refer to
    pbinom(c(13, 11, 10, 9, 7), size = 16, p = 9/11)
    pbinom(48, 62, 9/11)
    pbinom(9, 62, 1/11)
    pbinom(5, 62, 1/11)

    print('Finished Eight Runs (16 now, used to be 8!!)')

    ##
    ## Now produce outputs for all FFI run with the S_GCF model
    ##
    inds = grep('S_GCF', OrigSum[,1], fixed=TRUE)
    AllRunsCover = data.frame(Source=basename(OrigSum[inds,1]),
                              isEnveloped=OrigSum[inds,6],
                              UpBuf=OrigSum[inds, 7],
                              loBuf=OrigSum[inds,8])

    #####################################################
    ## OUTPUT FOR ALLEVENTS
    #ers = which(AllRunsCover[,2]==0)
    #FFI_no = unlist(lapply(strsplit(as.character(AllRunsCover$Source[ers]), '_'), f<-function(x) x[5]))
    #Mw_these = sapply(as.character(AllRunsCover$Source[ers]),f<-function(x) gsub('M',"", strsplit(x,'_')[[1]][3]) )
    #X = AllRunsCover[ers,]
    #allRuns_Err = (-X[,3]*(X[,3]<0) + X[,4]*(X[,4]<0))*10
    #
    #outTable = data.frame(FFI=FFI_no, Mw=as.character(Mw_these), DistanceErr=allRuns_Err)
    #write.table(outTable, file='AllRunsStats.txt', sep=" & ",
    #            row.names=FALSE, col.names=TRUE, quote=F,
    #            eol=' \\\\\n')


    # OUTPUT FOR ALL EVENTS, not just the ones with errors
    FFI_no = unlist(lapply(strsplit(as.character(AllRunsCover$Source), '_'), f<-function(x) x[5]))
    Mw_these = sapply(as.character(AllRunsCover$Source), f<-function(x) gsub('M',"", strsplit(x,'_')[[1]][3]) )
    X = AllRunsCover
    allRuns_Err = (-X[,3]*(X[,3]<0) + X[,4]*(X[,4]<0))*10
    reformatInds = (allRuns_Err!=0)
    allRuns_Err[reformatInds] = paste('\\textbf{', allRuns_Err[reformatInds], '}', sep="")
    
    allRuns_Err=paste(allRuns_Err, ' (',Orig_max[inds]*10, ')')

    outTable = data.frame(FFI=FFI_no, Mw=as.character(Mw_these), 
        DistanceErr=allRuns_Err)
    sort_Mw = sort(Mw_these,index.return=T)$ix
    outTable = outTable[sort_Mw,]

    newOutTable = cbind(outTable[1:21,], outTable[22:42,], 
        rbind(outTable[43:62,], data.frame(FFI="",Mw="",DistanceErr="")))
    write.table(newOutTable,file='AllRunsStats_includingSuccess.txt',sep=" & ",
                row.names=FALSE,col.names=TRUE,quote=F,
                eol=' \\\\\n')


    ####################################################
    maxSim = c()
    minSim = c()
    maxFFI = c()
    minFFI = c()
    nameOfSource = c()
    Events = grep('S_GCF', names(allWDL))
    for(i in 1:length(Events)){
        tmp = is_enveloped(allWDL[[Events[i]]], returnMinMax=T)
        maxSim = c(maxSim, tmp[1]) 
        minSim = c(minSim, tmp[2]) 
        maxFFI = c(maxFFI, tmp[3]) 
        minFFI = c(minFFI, tmp[4]) 
        nameOfSource = c(nameOfSource, names(allWDL)[Events[i]])
    }

    Mw = unlist(lapply(strsplit(basename(nameOfSource),'_'), 
        f<-function(x) as.numeric(gsub('M',"", x[[3]]))))

    pdf('InundationRatio.pdf', width=10, height=7)
    plot(Mw, maxSim/maxFFI, ylim=c(0.2,5.0), log='y', pch=15,
        xlab="", ylab="", cex.axis=1.3, cex=1.5)
    title(xlab=bquote(M[w]), cex.lab=2)
    title(ylab='Inundation Ratio',cex.lab=1.5,line=2.3)
    abline(h=1)
    points(Mw, minSim/minFFI,col=3,pch=19)
    legend('bottomright', 
        c('Max SFFM inundation / Max FFI inundation',
            'Min SFFM inundation / Min FFI inundation'),
        col=c(1,3), pch=c(15,19), cex=1.3, pt.cex=c(1.5,1))
    dev.off()


    # Export the summary table
    new_sumAWDL = sumAWDL
    new_sumAWDL$meanInun = (500-new_sumAWDL$meanInun)*10
    # NOTE: In the next 2 lines the names are the wrong way around (max=min,
    # min=max), and then they are fixed
    new_sumAWDL$maxInun = (500-new_sumAWDL$maxInun)*10
    new_sumAWDL$minInun = (500-new_sumAWDL$minInun)*10
    # Fix the names (min-->max, max--> min)
    myN = c(grep('maxInun', names(new_sumAWDL)), 
        grep('minInun', names(new_sumAWDL)))
    names(new_sumAWDL)[myN] = c('minInundation','maxInundation')

    write.table(new_sumAWDL, file='AllScenariosSummaryTable.csv', 
        sep=",", row.names=F)

    return(environment())
}

plotEnvelopeExamples<-function(k2Model, eventflag, titleFlag=1){
    # Make a plot of the real+synthetic tsunami inundation
    #
    
    if(!exists('allWDL')) allWDL=readRDS('all_wet_dry_lines.RDS')

    par(family='serif')
    eventInds=unlist(lapply( strsplit(basename(names(allWDL)), '_'), f<-function(x) as.numeric(x[5])))
    modelType=basename(dirname(names(allWDL)))
     
    for(k2 in k2Model){
        for(event in eventflag){
            eventI=which(grepl(k2, modelType) & eventInds==event)
            if(length(eventI)==0){
                print(paste0('No event ', event))
                next
            }
            allLines=allWDL[[eventI]]
            envelope=is_enveloped(allLines, returnSynthEnvelope=TRUE)        
            isEnv=is_enveloped(allLines)[1]
            #
            # Make plots
            #
            XX=envelope$x*10 # x coordinate in m
            # Enveloping polygon, in m
            envPoly=cbind(c(XX, rev(XX)), 
                          c((500-envelope$uLim)*10, (500-rev(envelope$lLim))*10))
            # FFI simulation in m
            origL=(500-envelope$orig)*10.
            #plot(range(XX), c(0, max(c(envPoly[,2], origL, 200))), 
            #     col=0,xlab="",ylab="")#,xlim=c(7000,9000))
            plot(c(0, max(c(envPoly[,2], origL, 200))), range(XX), 
                 col=0,xlab="",ylab="",cex.axis=2, yaxs='i')#,xlim=c(7000,9000))

            #title(xlab='Along-shore distance (m)')
            #title(ylab='Inundation distance (m)')
            title(ylab='Along-shore distance (m)',cex.lab=2,line=2.5)
            title(xlab='Inundation distance (m)',cex.lab=2)
            
            #polygon(envPoly,
            polygon(envPoly[,2:1],
                    col=rgb(200,200,200,alpha=140,maxColorValue=255),
                    border=NA)
            #points(envPoly,t='l',lwd=3,col='grey')
            points(envPoly[,2:1],t='l',lwd=3,col='grey')
            if(titleFlag==0){
                titleWords=paste0(basename(names(allWDL)[eventI]),' : ', as.character(isEnv))
                title(main=titleWords)
            }else{
                #title(main=titleFlag)
                Mw=gsub("M", "", strsplit(basename(names(allWDL)[eventI]), '_')[[1]][3])
                titleWords=bquote('FFI '~ .(event) ~ ' ('~M[w] ~'='~ .(Mw) ~')')
                title(main=titleWords,cex.main=3)
            }

            for(i in 1:length(allLines)){
                #points(allLines[[i]][,2]*10,(500-allLines[[i]][,1])*10,
                points((500-allLines[[i]][,1])*10,allLines[[i]][,2]*10,
                        t='l',lty='longdash',col='blue',lwd=0.5)
            }
            #points(XX,origL,t='l',lwd=2)
            points(origL,XX, t='l',lwd=3)

            abline(v=0,lty='dashed',col='red')
        }
    }

}

#makeInundationPdfs<-function(){
#    for (model in c('TSUNAMI_ABS','TSUNAMI_CLIP', 'TSUNAMI_ABS_TAPER','TSUNAMI_CLIP_TAPER')){
#        pdfname=paste0(model,'.pdf')
#        pdf(pdfname,width=5,height=10)
#        #pngname=paste0(model,'.png')
#        #png(pngname,width=5,height=10,res=300,units='in')
#        plotEnvelopeExamples(model, 1:66)
#        dev.off()
#
#    }
#}


makeNicePlot<-function(){

    #pdf('InundationFig.pdf', width=14,height=10)
    png('InundationFig.png', width=14,height=10,res=300,units='in')
    par(mfrow=c(1,4))
    par(mar=c(4,4,3,0.5))
    plotEnvelopeExamples('S_GCF', c(12, 44, 52, 9) )
    dev.off()

}


furtherInundationAnalysis<-function(){

    sum_AWDL = read.csv('AllScenariosSummaryTable.csv')

    # Create paths of okada displacement corresponding to each routine
    xyz_paths = lapply( strsplit(as.character(sum_AWDL[,2]), '/'), 
        f<-function(x) paste0(x[3], '/', x[4], '/', x[5], '.xyz'))
    xyz_paths = paste('../SLIP_MODELLING/OUTPUTS_2/DEFORMATION/', unlist(xyz_paths), sep="")
   
    # Find the peak okada displacement for each 
    range_okada = list()
    for(i in 1:length(xyz_paths)){
        range_okada[[xyz_paths[i]]] = range(matrix(scan(xyz_paths[i], sep=","), ncol=3, byrow=T)[,3])
    }
    peak_okada = unlist(lapply(range_okada, diff))

    # Categorize by model type 
    FFI_events = grepl('_Orig.xyz', xyz_paths)
    sffm_models = c('S_GA_', 'S_GAF', 'S_GC_', 'S_GCF', 
                    'S_SA_', 'S_SAF', 'S_SC_', 'S_SCF')

    model_type = rep(NA, length(sum_AWDL[,1]))
    for(i in 1:length(sffm_models)){
        model_type[grepl(sffm_models[i], xyz_paths)] = i
    }

    # Plot
    par(mfrow=c(3,3))
    plot(peak_okada, sum_AWDL$meanInun, log = 'xy', 
        col=FFI_events + 1, pch=model_type)

    for(i in 1:length(sffm_models)){
        inds = which(model_type==i)
        plot(peak_okada[inds], sum_AWDL$meanInun[inds], log='xy',
             col=FFI_events[inds] + 1, pch=model_type[inds])
    }

    event_name = unlist(lapply(strsplit(xyz_paths, '/'), f<-function(x) x[6]))
    event_number = unlist(lapply(strsplit(event_name, '_'), f<-function(x) x[5]))

    cors = list()
    for(en in unique(event_name)){
        keep = which(event_name==en)
        cors[[en]] = cor(sum_AWDL$meanInun[keep], peak_okada[keep], method='s')
    }

    Mw = unlist(lapply(strsplit(names(cors), '_'), f<-function(x) as.numeric(gsub('M', '', x[3])) ) )

    plot(Mw, unlist(cors))

    # Boxplots of mean inundation + peak okada range
    # Highlight the FFI in both cases
    inds = which(model_type==8)
    inds_FFI = which(model_type==8 & FFI_events)
    boxplot(sum_AWDL$meanInun[inds] ~ event_number[inds], log='y')
    boxplot(sum_AWDL$meanInun[inds_FFI] ~ event_number[inds_FFI], add=T, col='red', border='red')
    boxplot(peak_okada[inds] ~ event_number[inds], log='y')
    boxplot(peak_okada[inds_FFI] ~ event_number[inds_FFI], add=T, col='red', border='red')

    # Check for relations between the 'rank' of FFI mean inundation, compared
    # with the corresponding SFFM. Is this a function of Mw?
    FFI_rank = rep(NA, length(Mw))
    i = 0
    Mw = Mw*NA
    for(en in unique(event_name)){
        i = i+1
        Mw[i] = as.numeric(gsub("M", "", strsplit(en, '_')[[1]][3]))
        inds = which((model_type==4) & (event_name == en) & (!FFI_events))
        inds_FFI = which( (model_type==4) & (FFI_events) & (event_name == en))
        FFI_rank[i] = sum(sum_AWDL$meanInun[inds] < sum_AWDL$meanInun[inds_FFI])
    }

    plot(Mw, FFI_rank, main='Is there a relationship between FFI mean inundation rank and Mw?')
    FFI_rank_vs_Mw = lm(FFI_rank ~ Mw)
    summary(FFI_rank_vs_Mw)
    cor.test(FFI_rank, Mw, method='s')
    cor.test(FFI_rank, Mw, method='p')

    t.test(FFI_rank)
 
}
