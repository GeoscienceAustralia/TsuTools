#
# Extract events of interest from the SRCMOD webpages
#
# Event choice based on earlier analysis (NH), filenames from that work used
# as identification (see initial_event_list.txt)
#
# These are identified in the SRCMOD data and converted to useful formats
#


# Get metadata + raw web files
collateData=new.env()
source('collate_data.R',local=collateData)

# Find the date+Mw of our events 
our_event_list = scan('initial_event_list.txt', what=character())
theDate = lapply(as.list(our_event_list), 
                 myfun<-function(x) strsplit(x,'_')[[1]][1:3])
theDate = matrix(unlist(theDate), ncol=3, byrow=T)
year = as.numeric(theDate[,1])
month = match(theDate[,2], month.name)
day = as.numeric(theDate[,3])

dateCode = year*1e+06 + month*1e+03 + day
MwCode = unlist(
    lapply(as.list(our_event_list), 
           myfun<-function(x) as.numeric(strsplit(gsub('M',"",x),'_')[[1]][4])
          )
)

# NOTE: There is a double-up date code (2 Indonesia events on 2007-09-12). 
#       These can be distinguished by their Mw

# Convert the date from SRCMOD into a similar format
Sdates = collateData$metaAll$Date
Syear = unlist(lapply(as.list(Sdates), 
                      myfun<-function(x) as.numeric(strsplit(x,' ')[[1]][3])))
Sday = unlist(lapply(as.list(Sdates), 
                     myfun<-function(x) as.numeric(strsplit(gsub(',','',x),' ')[[1]][2])))
Smon1 = unlist(lapply(as.list(Sdates), 
                      myfun<-function(x) strsplit(x,' ')[[1]][1]))
Smon = match(Smon1,month.abb)
SdateCode = Syear*1e+06 + Smon*1e+03 + Sday
SMwCode = collateData$metaAll$Mw

# Get srsmod events which occur on the same day as our events
SRCMOD_MATCHES = which(SdateCode %in% dateCode)
#
# NOTE: There is a double-up date code (2 Indonesia events on 2007-09-12)
#       Manually deal with this now
SRCMOD_2_our_event_list = rep(NA,length(SRCMOD_MATCHES))
for(i in 1:length(SRCMOD_MATCHES)){
    # Hack to match the date, and if multiple dates match, then match the one
    # with the closest Mw
    matchCode = abs(SdateCode[SRCMOD_MATCHES[i]] + 
                    SMwCode[SRCMOD_MATCHES[i]]*0.1 - 
                    (dateCode+0.1*MwCode))
    SRCMOD_2_our_event_list[i] = which.min(matchCode)
}

# Check whether we could make rasters for all events matched. Sometimes we
# can't (e.g. due to multi-segment rupture models, etc) 
source('fsp2rast.R')
SRCMOD_rasts = plot_all(collateData$all_model_data)
# Search for events which match the date of Nick's events, but are NA (i.e.
# failed to parse)
SRCMOD_parse_failed = which(is.na(SRCMOD_rasts$rast_all[SRCMOD_MATCHES]))
# There are repeats in the SRCMOD rasters -- ensure we don't miss any events
missing_events = setdiff(1:length(our_event_list), 
                        SRCMOD_2_our_event_list[-SRCMOD_parse_failed])
stopifnot(length(missing_events)==0)

# Remove events for which we couldn't make a raster -- often these are 
# multi-segment ruptures [our code isn't designed for these], also there
# are some cases which look like data errors - see printed try-catch errors
SRCMOD_MATCHES = SRCMOD_MATCHES[-SRCMOD_parse_failed]
SRCMOD_2_our_event_list = SRCMOD_2_our_event_list[-SRCMOD_parse_failed]

# Remove events for which any dimensions are <= 5 cells
tooCoarse = which(unlist(lapply(SRCMOD_rasts$rast_all[SRCMOD_MATCHES], 
                              myfun<-function(x) min(dim(x)[1:2])<=5)))
SRCMOD_MATCHES = SRCMOD_MATCHES[-tooCoarse]
SRCMOD_2_our_event_list = SRCMOD_2_our_event_list[-tooCoarse]

# Get rasters which we could parse and were not removed above
SRCMOD_focus_rasts = plot_all(collateData$all_model_data[SRCMOD_MATCHES], 
                              'subduction_srcmod.pdf')

# Make some output files
outdir = 'SRCMOD_SUBDUCTION_DATA'
dir.create(outdir, showWarnings=FALSE)

write.table(collateData$metaAll[SRCMOD_MATCHES,], 
            paste0(outdir, '/SrcMod_MetaData.csv'), 
            row.names=FALSE, sep=",")

## Here we output a reduced srcmod table with the event tag
srcmod_tag_info=data.frame(Location=character(), 
                           Date=character(), 
                           Mw=numeric(), 
                           srcmodTag=character(), 
                           mytag=numeric())
for(i in 1:length(SRCMOD_MATCHES)){
    SMI = SRCMOD_MATCHES[i]
    eventi = collateData$metaAll$Location[SMI]
    datei = collateData$metaAll$Date[SMI]
    Mwi = collateData$metaAll$Mw[SMI]
    mytagi = i
    
    tagIndex = grep('EventTAG', collateData$all_model_data[[SMI]][[2]])
    srcmodTagi = strsplit(collateData$all_model_data[[SMI]][[2]][tagIndex], 
                          ' ')[[1]][3]

    srcmod_tag_info=rbind(srcmod_tag_info,
                          data.frame(Location=eventi, 
                                     Data=datei, 
                                     Mw=Mwi, 
                                     scrmodTag=srcmodTagi,
                                     mytag=mytagi))
}
write.table(srcmod_tag_info,
            file=paste0(outdir,'/SrcModTable_withTags.csv'),
            sep=",",row.names=FALSE)

## Below we output rasters
rastdir = paste0(outdir,'/slipRasts')
dir.create(rastdir,showWarnings=FALSE)

#for(i in 1:length(SRCMOD_focus_rasts$rast_all)){
for(i in 1:length(SRCMOD_MATCHES)){

    # Location name, simplified
    nameSimp = collateData$metaAll$Location[SRCMOD_MATCHES[i]]
    nameSimp = gsub(' ','-', nameSimp)
    nameSimp = gsub('(','', nameSimp,fixed=T)
    nameSimp = gsub(')','', nameSimp)
    nameSimp = gsub(',','-', nameSimp)
    # Author name, simplified
    authSimp = collateData$metaAll$Author[SRCMOD_MATCHES[i]]
    authSimp = gsub(' ','-',authSimp)
    authSimp = gsub('(','',authSimp,fixed=T)
    authSimp = gsub(')','',authSimp,fixed=T)
    authSimp = gsub(',','-',authSimp,fixed=T)
    authSimp = gsub('.','-',authSimp,fixed=T)
    
    # Name for raster -- make it informative + unique
    thisRastName = paste0(rastdir,'/','S_', SdateCode[SRCMOD_MATCHES[i]],'_M',
                          SMwCode[SRCMOD_MATCHES[i]],'_',nameSimp,'_',i,'.tif')

    # fsp2rast outputs a format slightly different to what we have used for the
    # random slip analysis.
    # Adjust the format here so that y increases down-dip, with the 'top' edge
    # having y=0
    m1 = SRCMOD_focus_rasts$rast_all[[i]]
    m1Mat = as.matrix(m1)
    m1MatFlip = m1Mat[(dim(m1)[1]):1,]
    m2 = raster(m1MatFlip,xmn=0, xmx=abs(extent(m1)@xmax-extent(m1)@xmin),
                ymn=0,ymx=abs(extent(m1)@ymax-extent(m1)@ymin))
    writeRaster(m2, filename=thisRastName,
                format='GTiff', options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

    # Name for FSP file
    thisFSPname = paste0(rastdir,'/','S_', SdateCode[SRCMOD_MATCHES[i]],'_M',
                         SMwCode[SRCMOD_MATCHES[i]],'_',nameSimp,'_',i,'.FSP')
    cat(collateData$all_model_data[[SRCMOD_MATCHES[i]]][[2]], 
        file=thisFSPname,sep="\n")
}

myrasts = Sys.glob(paste0(rastdir,'/*.tif'))
pdf('subduction_rastsFlipped.pdf')
for(myrast in myrasts){
    plot(raster(myrast))
    title(myrast)
}
dev.off()
