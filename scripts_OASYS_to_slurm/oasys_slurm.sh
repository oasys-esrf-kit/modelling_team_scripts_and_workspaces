#!/bin/bash -l
# Here we defined the python that we would like to use  
CMD_LAT=$@

CMD=$(echo /cvmfs/hpc.esrf.fr/software/packages/ubuntu20.04/x86_64/oasys/2023.11.23/bin/python \
       $CMD_LAT )

echo $CMD
$CMD
