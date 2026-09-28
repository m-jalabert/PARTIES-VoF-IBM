#!/bin/bash -l
#SBATCH --job-name=B14disc
#SBATCH --output=disc.out
#SBATCH --error=disc.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=01:00:00
# B.1.4 discriminator: identical to StageB_release2d but with ZERO particles.
# If the 9.3% meltwater conservation drift persists, it is in the CH/tracer time
# integration; if it vanishes, it is IBM-coupled.
set -euo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
C="${SLURM_SUBMIT_DIR}"
# Run output goes to SCRATCH, never $HOME (25 GB quota).  See roadmap B.0.1.
d="/anvil/scratch/x-mjalabert/ECCO_StageB/StageB_tracer_noparticle"
rm -rf "$d"; mkdir -p "$d"
cp "$C/parties.inp" "$C/stop.inp" "$C/p_fixed.inp" "$C/p_mobile.inp" "$d/"
cp /home/x-mjalabert/PARTIES/PARTIES/parties.stageB_rel5 "$d/parties"
cd "$d" && srun ./parties > run.log 2>&1
echo "tracer discriminator done"
