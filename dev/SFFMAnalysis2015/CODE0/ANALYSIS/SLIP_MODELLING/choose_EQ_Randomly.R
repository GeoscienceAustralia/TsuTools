# Choose a subset of the earthquakes 'randomly' but ensure they are all above
# 7.5, and have 4 (previously 2) in every 0.5 Mw increment above that.

chooseSources<-function(){
    # Choose FFI randomly, 16 in total, 8 at a time
    # The reason for this strange randomization method (i.e. 2 lots of 8 rather
    # than 1 lot of 16) is that initially I only ran 8 FFI -- but following
    # review + extension to more S* models, decided to extend to 16
    # models. To avoid changing the original set, we generate 2 sets of 8,
    # ensuring no overlap

    SourceNames = basename(Sys.glob('OUTPUTS_2/DEFORMATION/142769980501406_S_GA__NST_abs_SSD_none_RCS_TRUE_hurst_1_noise_gaussian/S_*'))
   
    Mw = unlist(lapply(strsplit(SourceNames, '_'), 
              f<-function(x) as.numeric(gsub('M',"",x[3]))))

    already_chosen_EQ_ind = c()
    set.seed(1) # Reproducible 
    for(repeats in 1:2){
        Mw_bins = seq(6.5,9.5,by=0.5)-0.0001
        Mw_bin = findInterval(Mw,Mw_bins)
        chosen_EQ_ind = c() #rep(NA,length(Mw_bins)-1)
        for(i in 1:max(Mw_bin)){
            # Skip events < 7.5
            if(i %in% c(1,2)) next 
        
            # Find events in the ith bin
            InBin = setdiff(which(Mw_bin==i), already_chosen_EQ_ind)

            # Select 2 at random
            chosen_EQ_ind = c(chosen_EQ_ind, sample(InBin,size=2))
        }

        write.table(SourceNames[chosen_EQ_ind], 
            file=paste0('subset_gt_7p5_', repeats,'.txt'), 
            row.names=F,col.names=F,quote=F)

        already_chosen_EQ_ind = chosen_EQ_ind
    }

}

copySources<-function(){
    # Code to copy the sources we need to a location, so we can put them on NCI
    # and RUN
    #
    # NOTE: The only reason this is coded in 2 sets of 8 events is that
    # initially only 1 set of 8 events was ran. To make it most simple to run
    # another 8, I chose to copy the 2nd set of 8 to a different directory

    source_subset_files = c('subset_gt_7p5_1.txt', 'subset_gt_7p5_2.txt')
    source_dest_folders = c('DEFORMATION_COPY', 'DEFORMATION_COPY2')

    for(ii in 1:length(source_subset_files)){

        chosenSources = scan(source_subset_files[ii], what='character')

        allSourceTypes = Sys.glob('OUTPUTS_2/DEFORMATION/*')

        for(src in chosenSources){
            for(srcType in allSourceTypes){
                allXYZ = Sys.glob(paste0(srcType,'/',src, '/*.xyz'))
                allXYZ_new = gsub('DEFORMATION', source_dest_folders[ii], allXYZ)
                dir.create(dirname(allXYZ_new[1]),recursive=TRUE, showWarnings=FALSE)
                for(i in 1:length(allXYZ)){
                    if(!file.exists(allXYZ_new[i])) file.copy(allXYZ[i], allXYZ_new[i])
                }
            }
        }
    }
}


copySingleSource<-function(){
    # Use to copy all the S_GCF sources which I have not previously put in
    # 'DEFORMATION_COPY and 'DEFORMATION_COPY2'
    ruptures_SGCF = Sys.glob('OUTPUTS_2/DEFORMATION/*S_GCF*/S_*')
    ruptures_already_done = c( Sys.glob('OUTPUTS_2/DEFORMATION_COPY/*S_GCF*/S_*'), 
        Sys.glob('OUTPUTS_2/DEFORMATION_COPY2/*S_GCF*/S_*'))

    dir.create('OUTPUTS_2/S_GCF', showWarnings=FALSE)

    ruptures_to_do = setdiff(basename(ruptures_SGCF), basename(ruptures_already_done))

    rupture_src_dir = dirname(ruptures_SGCF[1])

    for(i in 1:length(ruptures_to_do)){
        src_file = paste0(rupture_src_dir, '/', ruptures_to_do[i])
        dest_file = paste0('OUTPUTS_2/S_GCF/')
        print(src_file)
        file.copy(src_file, dest_file, recursive=TRUE)
    }
}
