# Read all the FFM data from Mai SRCMOD website
basedir='http://equake-rc.info/SRCMOD/searchmodels/allevents/'
extensiondirs=c('', paste0('?page=', 2:8))

# Read the base link page(s)
base_link_page=c()
for(i in 1:length(extensiondirs)){
    base_link_page=c(base_link_page,readLines(paste0(basedir,extensiondirs[i])))
}

# Extract the 'core' name for each FFM
model_link_ind=grep('SRCMOD/searchmodels/viewmodel', base_link_page)
link_splitter=strsplit(base_link_page[model_link_ind], '/')
model_corename=unlist(lapply(link_splitter, myfun<-function(x) x[5]))

model_fullnames=paste0("http://equake-rc.info/SRCMOD/searchmodels/viewmodel/", model_corename, '/')

fsp_data_basename='http://equake-rc.info/media/srcmod/_fsp_files/'

library(XML)
parse_page<-function(model_fullname){
    # read page text
    pagetxt=readLines(model_fullname)
    # get tables
    pagetables=readHTMLTable(pagetxt)
    # Find fsp file location
    data_file_line_index=grep('_fsp_files', pagetxt)
    fspname=strsplit(strsplit(pagetxt[data_file_line_index], '/')[[1]][5], '\"')[[1]][1]
    # File location 
    fsp_file_location=paste0(fsp_data_basename,  fspname)
    fsp_file=readLines(fsp_file_location)    
    return(list(pagetables, fsp_file, pagetxt))
}

# Download all the data
all_model_data=list()
for(i in 1:length(model_fullnames)){
    print(paste('fault', i, '...'))
    all_model_data[[i]] = parse_page(model_fullnames[i])
}

saveRDS(all_model_data,file='all_model_data.RDS')
