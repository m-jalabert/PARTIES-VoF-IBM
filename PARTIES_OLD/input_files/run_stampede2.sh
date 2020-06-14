#!/bin/bash
#SBATCH -J PARTIES_default         # Job name
#SBATCH -o log.out       # Name of stdout output file
#SBATCH -e log.err       # Name of stderr error file
#SBATCH -p skx-normal          # Queue name
#SBATCH -N 1               # Total # of nodes (now required)
#SBATCH -n 48             # Total # of mpi tasks
#SBATCH -t 1:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A TG-CTS150053

# Other commands must follow all #SBATCH directives...

module reset
module list
pwd
date

# Launch MPI job...
ibrun ./parties
