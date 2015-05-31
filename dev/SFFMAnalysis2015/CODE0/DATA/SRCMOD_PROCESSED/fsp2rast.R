## Rountines to convert FSP files to a raster format, check for various errors,
## make plots, etc

library(raster)
binme<-function(gridvals, Dx){
    ##########################################################################
    #
    # Function to 'bin' x or y values which are nearly on a grid (or supposed to be!)
    #
    # 'gridvals' is a set of x OR y coordinates, which in-theory should consist
    # of a finite set of values which are spaced by exactly 'Dx' (as we expect
    # for coordinates on a grid) 
    #
    # In practice for our dataset, there can be violations to the
    # ideal behaviour of 'gridvals' -- sometimes substantial ones.
    #
    # The current function is used to re-bin the gridvals appropriately
    
    # Define desirable grid values
    central_grid_vals = seq(min(gridvals), max(gridvals), by=Dx)

    # Sometimes round-off or deviations of 'gridvals' from regularity can cause
    # the seq to stop nearly Dx from max(gridvals). Fix this if it is bad
    # enough to break our binning
    if(max(gridvals)-max(central_grid_vals) > Dx/2){
        central_grid_vals = seq(min(gridvals), max(gridvals), 
            len = length(central_grid_vals) + 1)
    }

    # Double check that the max value of central_grid_vals really is
    # 'sufficiently close' to max(gridvals)
    # At least it must go in the correct bin.
    if(abs(max(gridvals)-max(central_grid_vals)) >= Dx/2){
        browser() 
        stop('Error in binning, should never occur')
    }

    # Make bin boundaries
    boundary_and_central_grid_vals = seq(min(gridvals)-Dx/2, max(gridvals)+Dx/2 , len=length(central_grid_vals)*2+1)
    mybreaks = boundary_and_central_grid_vals[ seq(1, length(boundary_and_central_grid_vals), 2) ]

    mybin = findInterval(gridvals, mybreaks)

    binned_values = central_grid_vals[mybin]
    if(max(abs(binned_values - gridvals)) > Dx/2){
        browser()
        stop("Error in binning -- coordinates don't seem close enough to being on a grid")
    }

    return(central_grid_vals[mybin])
}

#####################################################

fsp2rast<-function(fspfile){
    # Function to get FFI raster from FSP format
    # Takes an fspfile, outputs a raster with x=along-strike,y=down-dip, all in
    # kilometres. 
    # Note -- this raster format is slightly edited for our analysis in
    # 'produce_rasters_and_table.R'

    # Read the source in that txt format, and fix names
    fspfile_notab = gsub('\t', '%', fspfile) # Some files miss a comment char with a tab -- put it back
    myfsp = read.table(text=fspfile_notab, comment.char='%')

    # Get header
    myheaderRow = grep('  LAT  ', fspfile_notab)
    if(length(myheaderRow)!=1) stop('Failed to get one header')
    myheader = scan(text=fspfile_notab[myheaderRow], what='character', quiet=T)
    myheader = myheader[-1]
    myheader = gsub('X==EW', 'X', myheader)
    myheader = gsub('Y==NS', 'Y', myheader)
    names(myfsp) = myheader
    if(max(myfsp$Z)==0) stop('Z is all zero') # This happens sometimes

    if(min(diff(myfsp$Z))) stop('Z decreases -- probable data issue')


    # Get DX & DZ
    dxdzline = grep('Dx  =', fspfile) # Line with the data
    dx_dz_numeric = as.numeric(strsplit(fspfile[dxdzline], ' ')[[1]]) # Only non-NA values are Dx & Dz
    Dx = dx_dz_numeric[-which(is.na(dx_dz_numeric))][1]
    Dz = dx_dz_numeric[-which(is.na(dx_dz_numeric))][2]

    if(is.na(Dx)| is.na(Dz)){
        print('Problem getting Dx and/or Dz')
        print(paste('Dx : ', Dx, ' Dz ', Dz))
    }


    # Fit 2D plain to xyz data to check it is planar
    myplain = lm(myfsp$Z~myfsp$X+myfsp$Y)
    myplainsum = summary(myplain)
    if(myplainsum$r.squared<0.999){
        print('This might not be planar (or it has 90 deg dip)')
        print('Trying PCA solution -- careful')

        # PCA works well inthis case
        myplain = prcomp(cbind(myfsp$X, myfsp$Y, myfsp$Z))
        downcoord = myplain$x[,2]
        longcoord = myplain$x[,1]
    }else{
        # myplain defines a*x+b*y+constant=z    
        #
        # Coefficients of plain are related to normal vectors along-strike/down dip
        # As long as a*x+b*y = FIXED, then we are purely moving along strike
        # So, -b*x+a*y is a measure of position along strike (since it is
        # perpendicular to the line that is constant along strike). The actual
        # distances can be in km so long as (a,b) is normalised to have length 1.
        # 
        # Conversely, as long as -b*x+a*y = FIXED, we are purely moving down-dip. 
        # So, a*x+b*y is a measure of position down-dip. The distances would be km
        # AT THE SURFACE if (a,b) had length 1. But actually, we want distances
        # DOWN THE DIP. So we need (a,b,1) to have length 1. So we need to scale
        # by sqrt(a^2+b^2+1)/sqrt(a^2+b^2)
        a = coef(myplain)[2]
        b = coef(myplain)[3]

        unit_scale = sqrt(sum(c(a,b)^2))
        downcoord = (a*myfsp$X+b*myfsp$Y)/unit_scale*sqrt(a^2+b^2+1)
        longcoord = (-b*myfsp$X+a*myfsp$Y)/unit_scale
    }

    # Check that downcoord really does go down!
    if(cor(downcoord, myfsp$Z)<0.99){
        # Coordinates may be reversed
        if(cor(downcoord, myfsp$Z)< -0.99){
            downcoord = -downcoord
        }else{
            stop('POOR CORRELATION')
        }
    }

    # Normalise to downward coordinate to be '0' at the top
    downcoord = downcoord-max(downcoord)

    # Correct for round-off or non-planar behaviour with binning
    downcoordb = binme(downcoord,Dz)
    longcoordb = binme(longcoord,Dx)
    
    # Find unique values
    uniqueDown = sort(unique(downcoordb))
    uniqueLon = sort(unique(longcoordb))

    ##  -- check for bad behaviour (i.e. not regular grid)
    downward_spacing_range = range(diff(uniqueDown))
    along_spacing_range = range((diff(uniqueLon)))
    if(diff(downward_spacing_range) > 0.05*max(abs(downward_spacing_range))){
        print(uniqueDown)
        print(downward_spacing_range)
        stop('Irregular downward spacing')
    }
    if(diff(along_spacing_range) > 0.05*max(abs(along_spacing_range))){
        print(uniqueLon)
        print(along_spacing_range)
        stop('Irregular along-strike spacing')
    }

    # Associate each x/y coordinate with an index
    lonIndex = match(longcoordb,uniqueLon)
    downIndex = match(downcoordb,uniqueDown)    

    # Make a matrix
    mymat = matrix(NA, ncol=length(uniqueLon), nrow=length(uniqueDown))
    mymat[cbind(downIndex,lonIndex)] = myfsp$SLIP

    if(prod(dim(mymat))!=length(downcoord)){
        stop('Different number of raster cells and points')
    }

    myxcell = mean(diff(uniqueLon))
    mydowncell = mean(diff(uniqueDown))
   
    myrast=raster(mymat,
                  xmn=min(uniqueLon)-0.5*myxcell,
                  xmx=max(uniqueLon)+0.5*myxcell,
                  ymn=min(uniqueDown)-0.5*mydowncell,
                  ymx=max(uniqueDown)+0.5*mydowncell)

    return(myrast)
}

##############################################
plot_all<-function(all_model_data, pdfname='quickplots.pdf'){
    ## Plot all the rasters (if they don't have errors), and 
    ## save them in a list in the function environment (which is returned)
    pdf(pdfname, width=8, height=6)
    count_errs = 0
    rast_all = list()
    for(i in 1:length(all_model_data)){
        print(i)
        m1 = try(fsp2rast(all_model_data[[i]][[2]]))
        if(class(m1)=='try-error'){
            count_errs = count_errs+1
            rast_all[[i]]=NA
            next
        }
        plot(m1,asp=1)
        mystring = paste0(as.character(all_model_data[[i]][[1]][[1]][[2]]), ' ',
                          as.character(all_model_data[[i]][[1]][[1]][[3]]), ' ',
                          as.character(all_model_data[[i]][[1]][[2]][[1]][2]))
        title( mystring)
        rast_all[[i]] = m1
    }
    dev.off()
    return(environment())
}

