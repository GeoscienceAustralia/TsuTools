Before these routines are run, the data must be pre-processed with the code in ../DATA/SRCMOD_PROCESSED

Then the codes here can be run by starting R in this directory, and doing:

    source('driver.R')

Package dependencies must be installed beforehand. Non-base-R packages include
    - raster, rgdal, igraph, HyperboicDist, stabledist, fitdistrplus [on CRAN]
    - EqSim, colorRampPC [in the 'packages' directory]

Other codes of relevance [not called by driver.R]:

    calc_kcx_kcy_1D.R implements some 1D fitting methods for the corner wavenumbers (for comparison)

    check_hurst.R investigates fitting models with variable hurst exponent
    (from which we concluded the fit is too variable to be worth proceeding with)

    choose_EQ_Randomly.R was used to generate various random FFI subsets and take care of copy operations.
    
    stable_distribution_fits_checks.R fits a stable distribution to the 'white-noise'
    of each FFI + contains code to check the fitting method. This method was not
    used for the final paper because it does not account for the changes that the SFFM
    causes to the phase, and also the so called 'white-noise' we extract from the SFFM
    is not iid. However, it contains some useful tests confirming our stable distribution
    fitting routines.

    more_stable_alpha_checks.R includes better code to estimate the
    stable-distribution alpha parameter for SFFM (as described in the
    Supplementary material for the paper) 

    


