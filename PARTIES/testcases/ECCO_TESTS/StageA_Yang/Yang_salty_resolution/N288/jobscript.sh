#!/bin/bash -l
#SBATCH --job-name=YangSalt288
#SBATCH --output=salt288.out
#SBATCH --error=salt288.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=6:00:00
#SBATCH --mail-user=mjalabert@ucsb.edu
#SBATCH --mail-type=END,FAIL

# Salty Sm=5 dSv=5 at N=288, t -> 350.0.  Targets: Yang ~216, PARTIES 1440^2 = 179.85.

set -euo pipefail
module load gcc/11.2.0 openmpi/4.1.6 fftw/3.3.8

PARTIES_BIN="${PARTIES_BIN:-$HOME/PARTIES/PARTIES/parties}"
if [ ! -x ./parties ]; then cp "$PARTIES_BIN" ./parties; fi

echo "Yang salty N=288 on $SLURM_NTASKS ranks, started $(date)"
srun ./parties > run.log 2>&1
echo "finished $(date)"
