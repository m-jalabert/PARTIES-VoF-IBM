#!/bin/bash -l
#SBATCH --job-name=SBnoninvF
#SBATCH --output=noninv_fresh.out
#SBATCH --error=noninv_fresh.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=01:00:00
# NOTE: with LAG_PARTICLE_RESOLVED enabled PARTIES requires p_fixed.inp and
# p_mobile.inp in the run directory even when there are zero particles
# (ParticleInput.c:87).  Every Stage-B deck must ship them.
set -euo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
C="${SLURM_SUBMIT_DIR}"
d="$C/fresh288_stageB"; mkdir -p "$d"
cp "$C/parties_fresh288.inp" "$d/parties.inp"
cp "$C/stop.inp" "$C/p_fixed.inp" "$C/p_mobile.inp" "$d/"
cp /home/x-mjalabert/PARTIES/PARTIES/parties.stageB "$d/parties"
cd "$d" && srun ./parties > run.log 2>&1
echo "fresh288 done"
