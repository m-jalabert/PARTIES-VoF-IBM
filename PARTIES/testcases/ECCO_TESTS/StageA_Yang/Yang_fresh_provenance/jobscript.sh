#!/bin/bash -l
#SBATCH --job-name=YangProv1440
#SBATCH --output=prov1440.out
#SBATCH --error=prov1440.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=10:00:00
#SBATCH --mail-user=mjalabert@ucsb.edu
#SBATCH --mail-type=END,FAIL

# A.5.2 provenance re-run: 1440^2 Sm=0, t -> 10, current binary.
# 1440 = 2^5*3^2*5; 128 = 16x8 -> 90x180 cells per rank.
# CAVEAT: the completed run used 256 ranks.  A PASS therefore exonerates both
# the binary and the rank count at once; a FAIL would need one follow-up to
# separate them.  Rank-independence is a stated design requirement and the
# ladder already spans 64/128 ranks on a smooth trend, so this is low risk.

set -euo pipefail
module load gcc/11.2.0 openmpi/4.1.6 fftw/3.3.8

PARTIES_BIN="${PARTIES_BIN:-$HOME/PARTIES/PARTIES/parties}"
if [ ! -x ./parties ]; then cp "$PARTIES_BIN" ./parties; fi

echo "Yang 1440^2 provenance re-run on $SLURM_NTASKS ranks, started $(date)"
srun ./parties > run.log 2>&1
echo "finished $(date)"
