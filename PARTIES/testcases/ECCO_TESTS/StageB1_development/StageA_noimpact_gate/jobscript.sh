#!/bin/bash -l
#SBATCH --job-name=SAnoimpact
#SBATCH --output=gate.out
#SBATCH --error=gate.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:30:00
# Proves the ECCO additions do not touch the pre-ECCO code path: Stage-A flags
# OFF, run the recorded Le=100 salt gate, compare against lambda = 0.21844955.
# NOTE: Stage-A build has LAG_PARTICLE_RESOLVED off, so NO p_*.inp are needed.
set -euo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
C="${SLURM_SUBMIT_DIR}"
# Run output goes to SCRATCH, never $HOME (25 GB quota).  See roadmap B.0.1.
d="/anvil/scratch/x-mjalabert/ECCO_StageB/StageA_noimpact_gate"
rm -rf "$d"; mkdir -p "$d"
cp "$C/parties.inp" "$C/stop.inp" "$d/"
cp /home/x-mjalabert/PARTIES/PARTIES/parties.stageA_recheck "$d/parties"
cd "$d" && srun ./parties > run.log 2>&1
echo "stageA no-impact gate done"
