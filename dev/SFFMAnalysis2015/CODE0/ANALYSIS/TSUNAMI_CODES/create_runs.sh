#!/bin/bash
#
# Script to create ANUGA runs for all deformations, from template
#

#SLIPMODEL='142769980501406_S_GA__NST_abs_SSD_none_RCS_TRUE_hurst_1_noise_gaussian'

# NOTE: This script requires manual adjustment for different sets of runs -- to the
# directoryname where the slip models are held, and the name where the ANGUA models are run

for SLIPMODEL in S_GCF/*; 
    do SM=$(basename $SLIPMODEL);
    RUNDIR='BEST_MODEL_REMAINING_FFI/'$SM
    DEFORMATIONDIR='S_GCF/'$SM
    mkdir -p $RUNDIR
    for i in $DEFORMATIONDIR/S_*; 
        do bi=$(basename $i);
        mkdir -p $RUNDIR/$bi;
        echo $i;
        for j in $i/*.xyz;
            do bj=$(basename $j .xyz) 
            mkdir -p $RUNDIR/$bi/$bj;
            echo $j;
            # Copy codes
            cp -r TEMPLATE_CODES/* $RUNDIR/$bi/$bj
            # Copy initial condition
            cp $j $RUNDIR/$bi/$bj/INPUT_DATA/
            # Fix file reference in project.py
            sed -i s/REPLACEWITHSED/$bj'_SMOOTHED.xyz'/g $RUNDIR/$bi/$bj/project.py
            done
    done
done
