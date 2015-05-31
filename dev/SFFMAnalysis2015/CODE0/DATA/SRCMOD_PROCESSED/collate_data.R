## Script to extract lots of metadata (event site / geometry / etc) from SRCMOD
## webpage data (as stored in .RDS file)

all_model_data = readRDS('../SRCMOD_RAW/all_model_data.RDS')

get_event_data<-function(all_model_data){
    ###########################################################################
    # Extract FFM coordinates / dates/ authors /  names etc from all_model_data
    # store in 'metadata_event'
    metadata_event = all_model_data[[1]][[1]][[1]]
    for(i in 2:length(all_model_data)){
        metadata_event = rbind(metadata_event, all_model_data[[i]][[1]][[1]])
    }
    metadata_event[,1] = as.character(metadata_event[,1])
    metadata_event[,2] = as.character(metadata_event[,2])
    metadata_event[,3] = as.character(metadata_event[,3])

    for(i in c(4,5,6)){
        metadata_event[,i] = as.numeric(as.character(metadata_event[,i]))
    }
    metadata_event[,7] = as.character(metadata_event[,7])

    return(metadata_event)
}


fault_geom<-function(atext){
    # Hacky text parsing utility function for get_metadata_geometry
    # Get some fault geometric pars

    # Remove double tab
    l1 = gsub("\t\t", "\t", atext[7:9])
    l2 = gsub("\t", " ", l1)
    l3 = gsub("  ", " ", l2)
    l4 = unlist(strsplit(l3," "))

    mynumind = grep('=', l4)+1
    mynums = matrix(as.numeric(l4[mynumind]), nrow=1,byrow=T)
    colnames(mynums) = l4[mynumind-2]
    return(mynums)
}

get_metadata_geometry<-function(all_model_data){
    #####################################################################
    # Extract crude fault geometric parameters from all_model_data.
    # Store in 'metadata_geometry'
    
    metadata_geometry = fault_geom(all_model_data[[1]][[2]])
    namestore = names(metadata_geometry)
    for(i in 2:length(all_model_data)){
        metadata_geometry = rbind(metadata_geometry,fault_geom(all_model_data[[i]][[2]]))
    }
    metadata_geometry = as.data.frame(metadata_geometry)

    return(metadata_geometry)

}

# Combine metadata_event / metadata_geometry into a data.frame "metaAll"
metaAll = cbind(get_event_data(all_model_data),
                get_metadata_geometry(all_model_data))

