#!/bin/bash
#PBS -P n74
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -lmem=8GB
#PBS -lncpus=16
#PBS -l wd

## First batch of 8
#headDir='EIGHT_FFI/142769980501406_S_GA__NST_abs_SSD_none_RCS_TRUE_hurst_1_noise_gaussian'
#headDir='EIGHT_FFI/142779441194322_S_SA__NST_abs_SSD_none_RCS_TRUE_hurst_1_noise_stable'
#headDir='EIGHT_FFI/142771226860869_S_GAF_NST_abs_SSD_gaussian_RCS_TRUE_hurst_1_noise_gaussian'
#headDir='EIGHT_FFI/142772758403118_S_GC__NST_clip_SSD_none_RCS_TRUE_hurst_1_noise_gaussian'
#headDir='EIGHT_FFI/142778164514969_S_GCF_NST_clip_SSD_gaussian_RCS_TRUE_hurst_1_noise_gaussian/'
#headDir='EIGHT_FFI/142780735004539_S_SAF_NST_abs_SSD_gaussian_RCS_TRUE_hurst_1_noise_stable'
#headDir='EIGHT_FFI/142782217353858_S_SC__NST_clip_SSD_none_RCS_TRUE_hurst_1_noise_stable'
#headDir='EIGHT_FFI/142783390629079_S_SCF_NST_clip_SSD_gaussian_RCS_TRUE_hurst_1_noise_stable'

## 2nd batch of 8
#headDir='EIGHT_FFI2/142769980501406_S_GA__NST_abs_SSD_none_RCS_TRUE_hurst_1_noise_gaussian'
#headDir='EIGHT_FFI2/142779441194322_S_SA__NST_abs_SSD_none_RCS_TRUE_hurst_1_noise_stable'
#headDir='EIGHT_FFI2/142771226860869_S_GAF_NST_abs_SSD_gaussian_RCS_TRUE_hurst_1_noise_gaussian'
#headDir='EIGHT_FFI2/142772758403118_S_GC__NST_clip_SSD_none_RCS_TRUE_hurst_1_noise_gaussian'
#headDir='EIGHT_FFI2/142778164514969_S_GCF_NST_clip_SSD_gaussian_RCS_TRUE_hurst_1_noise_gaussian/'
#headDir='EIGHT_FFI2/142780735004539_S_SAF_NST_abs_SSD_gaussian_RCS_TRUE_hurst_1_noise_stable'
#headDir='EIGHT_FFI2/142782217353858_S_SC__NST_clip_SSD_none_RCS_TRUE_hurst_1_noise_stable'
#headDir='EIGHT_FFI2/142783390629079_S_SCF_NST_clip_SSD_gaussian_RCS_TRUE_hurst_1_noise_stable'

## Remaining models
headDir='BEST_MODEL_REMAINING_FFI/142778164514969_S_GCF_NST_clip_SSD_gaussian_RCS_TRUE_hurst_1_noise_gaussian'

maxSites=1

mydir=$(pwd)
echo $mydir
counter=0

for i in $(ls $headDir'/');
    do echo $i ;
    cd $mydir/$headDir/$i;
    if [ -a 'running_now.txt' ] 
        then 
          echo 'skipping ' $i; 
    else
	counter=$(($counter+1))    
        touch 'running_now.txt'
        for j in OceanInitial_*;
            do cd $mydir/$headDir/$i/$j;
        	echo $j ': running subThis.PBS'
        	qsub subThis.PBS
        	cd ..
        done
        # Break when we have done maxRuns
        if [ $counter -eq $maxSites ]
          then 
	     # Stop setting off runs
             cd $mydir
             break;	
        fi
    fi
    cd ../../
done

