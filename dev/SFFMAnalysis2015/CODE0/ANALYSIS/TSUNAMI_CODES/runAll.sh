#!/bin/bash
maxRuns=20


mydir=$(pwd)
counter=0
for i in RUNS/*/*/; 
    do cd $i
    if [ -a 'MODEL_OUTPUTS' ] 
        then 
          echo 'skipping ' $i; 
    else
        qsub subThis.PBS
	#echo 'would run this';
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
