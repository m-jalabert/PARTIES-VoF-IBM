#!/bin/bash -l
#SBATCH --job-name=B14fix
#SBATCH --output=fix.out
#SBATCH --error=fix.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=03:00:00
#
# Roadmap B.1.4 -- FIX arm.  Byte-identical deck to the control; the binary
# differs by exactly one compile flag, VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE, which
# returns liquid trimmed inside the resolved grain to the diffuse rim around it
# instead of discarding it.
#
# Pass criteria vs the control:
#   * SED_CLIP_AUDIT cum|sed_resid| stays ~0  (the rim absorbed the mass locally;
#     a large residual means it had to be teleported by the global band restore)
#   * the meltwater-tracer identity int C_mw = int dF closes to <= 1e-3
#     (control reproduces the recorded +9.3e-2 at t=30)
#   * release time and settling trajectory are not materially disturbed
set -euo pipefail
set +u; module purge; module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8; set -u
C="${SLURM_SUBMIT_DIR}"
# Run output goes to SCRATCH, never $HOME (25 GB quota).  See roadmap B.0.1.
d="/anvil/scratch/x-mjalabert/ECCO_StageB/B14_sedclip_fix_${SLURM_JOB_ID}"
mkdir -p "$d"
cp "$C/parties.inp" "$C/stop.inp" "$C/p_fixed.inp" "$C/p_mobile.inp" "$d/"
cp /home/x-mjalabert/PARTIES/PARTIES/parties.stageB_sedfix "$d/parties"
cd "$d" && srun ./parties > run.log 2>&1
echo "B.1.4 fix arm done"
