#!/bin/bash
maxRuns=100


mydir=$(pwd)
counter=0
for i in TSUNAMI_CLIP_TAPER/*/*/; 
    do cd $i
    # The following lines look for a merged sww file
    shopt -s nullglob
    M=(MODEL_OUTPUTS/*/earthquake*_earthquake.sww) 
    #echo $M
    #if [ -e 'MODEL_OUTPUTS' ]
    if [ "${#M[*]}" -ge 1 ]
        then 
	counter=$counter;
        #echo '.'
  	#echo 'skipping ' $i; 
    else
	echo 'Running ' $i;
        qsub subThis.PBS
	counter=$(($counter+1))    
        echo 'counter = ' $counter
        if [ $counter -eq $maxRuns ]
          then 
	     # Stop setting off runs
             cd $mydir
             break;	
        fi
    fi
    cd $mydir
    done
