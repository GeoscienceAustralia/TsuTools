##############################################################################

./DATA 
    contains scripts to download / preprocess key data -- see the readme there for more information

./PACKAGES 

    contains some of our R packages which are required for the analysis
    but not available on CRAN -- these need to be installed for the analysis,
    typically with the terminal command 'R CMD INSTALL desired_package_name' from
    inside the package directory.

    Several other packages are required too, but they are 'standard' CRAN packages
    which can be installed from within R with (e.g. to get the 'raster' package):
        install.packages('raster')

    If you run the code and it complains about not being able to find a package,
    you should be able to install it with a command like that above.
    
./ANALYSIS 
    contains the analyses

###############################################################################

To run the analysis, you would need to run the codes in the following
directories in order, as described in their README.txt files

1) ./DATA/SRCMOD_RAW and SRCMOD_PROCESSED 
    (to preprocess the data)
2) ./ANALYSIS/SLIP_MODELLING 
    ( main earthquake analysis )
3) ./ANALYSIS/TSUNAMI_CODES 
    ( To make the tsunami codes to be executed in a super-computing environment)
    -- Super computing resources are required to run all the tsunami simulations
    -- This involves some manual changes to paths in run-scripts.
    -- The following steps assume that is already done
4) ./ANALYSIS/TSUNAMI_ANALYSIS 
    ( To extract key information from the tsunami simulations )
5) ./ANALYSIS/PAPER_FIGURES 
    (To make the figures from the paper)
