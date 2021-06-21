#!/bin/bash

for kappa in 1875 #1800 1825 #1875
do
    folder=/home/pietro/Desktop/DATA_LOCAL/TEK_flow/nfh/n289b0350k5hf$kappa
    # folder=HYDRA/nfh/n169b0350k5hf$kappa
    

    for file in $(ls $folder)
    do
        jobn=$(echo $file | awk -F'_' '{print $2}' | sed 's/^0*//')
        if [ "$jobn" -lt "10" ]; then
        #     mv $file $folder/job_00$jobn
            mv $folder/$file $folder/job_00$jobn
        elif [ "$jobn" -lt "100" ]; then
            mv $folder/$file $folder/job_0$jobn
        fi
    done

    # for file in $(ls $folder/job_00*)
    # do
    #     jobn=$(echo $file | awk -F'_00' '{print $2}')
    #     mv $file $folder/job_$jobn
    # done

    
done