#!/bin/bash -l
#SBATCH --job-name=B13fn
#SBATCH --output=b13fn.out
#SBATCH --error=b13fn.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=04:00:00
# B.1.3 functional release test.  Every Stage-B deck must ship p_fixed.inp and
# p_mobile.inp (ParticleInput.c:87) and stop.inp (Temporal_int.c:495).
set -euo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
C="${SLURM_SUBMIT_DIR}"
# Run output goes to SCRATCH, never $HOME (25 GB quota).  See roadmap B.0.1.
d="/anvil/scratch/x-mjalabert/ECCO_StageB/StageB_release_fn"
rm -rf "$d"; mkdir -p "$d"
cp "$C/parties.inp" "$C/stop.inp" "$C/p_fixed.inp" "$C/p_mobile.inp" "$d/"
cp /home/x-mjalabert/PARTIES/PARTIES/parties.stageB_rel4 "$d/parties"
cd "$d" && srun ./parties > run.log 2>&1
echo "B13 functional done"
