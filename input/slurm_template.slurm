#!/bin/bash
     
#SBATCH --job-name=JOBNAME
#SBATCH --nodes=NODES
#SBATCH --ntasks=NTASKS
#SBATCH --cpus-per-task=NCPUS
#SBATCH --time=infinite
#SBATCH --nodelist nostromo-0-NODENUMBER

exit_handler() {
    date 
    echo "signal caught by shell script" 
    kill -SIGINT "${child}" 
}

trap exit_handler SIGINT INT TERM SIGTERM

 
# go to working directory
workdir=WORKDIR
cd $workdir

# Load necessary modules
module load foss


#run
../../../ReversibleSim_auto parameters.dat > log.out 2>err.out &
child=$!
wait "$child"
