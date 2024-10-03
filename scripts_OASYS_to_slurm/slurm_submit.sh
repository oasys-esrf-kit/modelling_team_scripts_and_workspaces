#!/bin/bash -l
#
# Usage ./slurm_submit.sh CONFIGURATION_FILE(S)
#

if [ $# -eq 0 ]
then
   echo "Error: no filename given"
   exit 1
fi

FILE=$@
echo $FILE
sbatch -p nice --nodes=1 --mincpus=28 --time=12:00:00 --job-name=$FILE \
   /home/esrf/user/etc/oasys_slurm.sh \
   $FILE
