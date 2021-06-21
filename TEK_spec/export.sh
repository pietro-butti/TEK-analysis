#!/bin/bash

destination=/home/pietro/Desktop/DATA_LOCAL/TEK_spec/$1
folder=/lustre/tekproy/pietro/mesons/GEVP/Analysis/$1

# mkdir DATA/$1
mkdir $destination

for op in 8 #4 5 6 7 8 9 10
do
    mpcac=$(echo $folder/Op$op'_t22_t11/*mpcac*dl')
    vi=$(echo $folder/Op$op'_t22_t11/*vi*cxdl')
    pp=$(echo $folder/Op$op'_t22_t11/*pp*cxdl')


    # destination=DATA/$1/Op$op'_t22_t11'
    destination_op=$destination/Op$op'_t22_t11'
    mkdir $destination_op

    scp tekproy@hydra.ift.uam-csic.es:$mpcac $destination_op
    scp tekproy@hydra.ift.uam-csic.es:$vi $destination_op
    scp tekproy@hydra.ift.uam-csic.es:$pp $destination_op

done



# T5kPr0y_2011?