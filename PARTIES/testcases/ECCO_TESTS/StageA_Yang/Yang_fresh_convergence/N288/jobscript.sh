#!/bin/bash -l
#SBATCH --job-name=YangFresh288
#SBATCH --output=fresh288.out
#SBATCH --error=fresh288.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=3:00:00
#SBATCH --mail-user=mjalabert@ucsb.edu
#SBATCH --mail-type=END,FAIL

# Roadmap A.5.2 freshwater convergence ladder, N = 288, t -> 140.0.
# Compare t_half against the completed 1440^2 point (60.287) and the digitized
# Yang 2-D dSv=0 reference (107.42).

set -euo pipefail

# Anvil toolchain.  openmpi/4.1.6 is what the binary actually links against
# (check with `ldd parties`); HDF5 and HYPRE come from ~/Software/PARTIES_Libs
# via the binary's RPATH, so they need no module.  Do NOT use the gnu14/openmpi5
# names -- they do not exist on this machine and `set -e` will kill the job at
# the module load.
module load gcc/11.2.0 openmpi/4.1.6 fftw/3.3.8

PARTIES_BIN="${PARTIES_BIN:-$HOME/PARTIES/PARTIES/parties}"
if [ ! -x ./parties ]; then
    cp "$PARTIES_BIN" ./parties
fi

echo "Yang freshwater N=288 on $SLURM_NTASKS ranks, started $(date)"
srun ./parties > run.log 2>&1
echo "finished $(date)"
