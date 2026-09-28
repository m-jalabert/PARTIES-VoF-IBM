#!/bin/bash -l
#SBATCH --job-name=Yang_production
#SBATCH --output=production_job.out
#SBATCH --error=production_job.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=36
#SBATCH --time=149:00:00
#SBATCH --mail-user=mjalabert@ucsb.edu
#SBATCH --mail-type=END,FAIL

# A.4 Yang production run: TWOD_CARTESIAN, 1440 x 1440 x 1 storage slab,
# published window t -> 200.
# Adjust nodes/ntasks to the cluster; total ranks should factor cleanly into
# an NPX x NPY x 1 grid (the code auto-factorizes; 1440 = 2^5*3^2*5, so e.g.
# 64, 100, 144, 225, 256, 400 ranks all decompose well with NPZ = 1).
#
# The run may not finish t = 200 in one wall-clock window. To continue:
# set "resume = 1" in parties.inp and resubmit -- the code restarts from the
# latest Resume.h5. Keep everything else unchanged.

set -euo pipefail

module load gnu14/14.2.0 hwloc/2.12.0 ucx/1.18.0 libfabric/1.18.0 openmpi5/5.0.7 hdf5/1.14.6

# Copy the (already built, full-Yang-flag) binary next to the inputs on first launch.
PARTIES_BIN="${PARTIES_BIN:-$HOME/PARTIES/PARTIES/parties}"
if [ ! -x ./parties ]; then
    cp "$PARTIES_BIN" ./parties
fi

echo "Starting Yang production on $SLURM_NTASKS ranks"
mpirun -np "$SLURM_NTASKS" ./parties > run.log 2>&1
echo "Yang production run completed (or hit wall clock)"
